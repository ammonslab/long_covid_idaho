# Author: Jeremy Boyd (jeremy.boyd@va.gov)
# Description: Simpler long covid model

# Organize data ---------------------------------------------------------------

# File listing sequences that failed Nextclade quality control
failed_qc <- read_csv("data/bad_qc.csv") |>
  mutate(
    seqName = factor(seqName),
    failed_qc = 1L
  ) |>
  select(seqName, failed_qc)

# CDW data
data <- read_feather("data/cdw_data.feather") |>

  # Better names for variables we're using
  rename(
    index_date = IndexDate,
    age = AgeAtIndexDate,
    sex = Sex,
    race = Race,
    ethnicity = Ethnicity,
    vax = n_vax_doses_prior,
    rurality = cdc,
    svi_num = RPL_THEMES,
    comorbidity = CCI2yrs,
    smoking = Smoke2yrs,
    paxlovid = paxlovid_flag,
    dex = dex_flag,
    no_ins = NoRecordOfInsurance
  ) |>
  mutate(
    # Patient is a factor
    patient = factor(patient),

    # Compute number of phenotypes
    n_pheno = neuro_flag + cardio_flag + dig_flag,

    # Compute outcome label
    outcome = case_when(
      n_pheno > 1 ~ "Multisystem",
      neuro_flag == 1 ~ "Neuropsychiatric",
      cardio_flag == 1 ~ "Cardiopulmonary",
      dig_flag == 1 ~ "Digestive",
      TRUE ~ "Recovered"
    ),
    outcome = fct_relevel(outcome, "Recovered"),

    # Convert NA race/ethnicity to "unknown"; standardize case
    across(c(race, ethnicity), ~ if_else(is.na(.x), "Unknown", .x)),
    across(c(race, ethnicity), ~ str_to_title(.x)),

    # Convert logicals to 1/0
    across(where(is.logical), ~ if_else(.x == TRUE, 1L, 0L)),

    # Convert age to double
    age = as.numeric(age),

    # Not finding any transformations that improve CCI normality, so I'm
    # instead defining a comorbidity flag for any case with CCI >= 1, which
    # is a common practice.
    comorbidity = factor(if_else(comorbidity >= 1, "Yes", "No")),

    # Change Smoke2yrs to be a flag indicating evidence of current/former
    # smoker status in the 2 years prior to index.
    smoking = factor(if_else(
      smoking %in% c("Current Smoker", "Former Smoker"),
      "Current/Former",
      "Non-Smoker"
    )),

    # Convert dates to dates
    across(c(index_date, matches("DateTime")), ~ as.Date(.)),

    # Compute flag indicating whether patient died within 270 days of index
    index_death_lag_days = as.numeric(difftime(
      DeathDateTime,
      index_date,
      units = "days"
    )),
    death_flag = if_else(index_death_lag_days <= 270, 1L, 0L),

    # Define clade groups
    clade = factor(case_when(
      clade_nextstrain %in%
        c("20A", "20B", "20C", "20G", "20I", "20J", "21C", "21H") ~ "Pre-Delta",
      clade_nextstrain %in% c("21A", "21I", "21J") ~ "Delta",
      clade_nextstrain %in% c("21M", "21K", "21L") ~ "Early Omicron",
      clade_nextstrain %in%
        c(
          "22A",
          "22B",
          "22C",
          "22D",
          "22E",
          "22F",
          "recombinant"
        ) ~ "Omicron 2022+",
      TRUE ~ clade_nextstrain
    )),

    # Define race
    race = factor(
      if_else(
        race %in%
          c(
            "Asian",
            "American Indian Or Alaska Native",
            "Native Hawaiian Or Other Pacific Islander",
            "Black Or African American"
          ),
        "Non-White",
        race
      )
    ),

    # Define rurality
    rurality = factor(if_else(rurality %in% c("1", "2", "3"), "Low", "High")),

    # Definition of SVI categories based on visual inspection of
    # numeric SVI distribution.
    svi = if_else(svi_num > 0.35, "High", "Low"),

    # Define vaccine
    vax = factor(if_else(vax %in% c("0", "1"), "< 2", ">= 2")),

    # Paxlovid labels
    paxlovid = factor(if_else(paxlovid == 1, "Yes", "No")),

    # Make any character predictors factor
    across(where(is.character), ~ factor(.x))
  )

# Compute number of missing values per row for select variables
missing <- data |>
  select(
    patient,
    index_date,
    # Outcomes
    outcome,
    # Demographics
    sex,
    age,
    comorbidity,
    smoking,
    # Treatment
    # NOTE: n_vax_does_prior is based on all Covid vaccine doses in the
    # patient's record, even going back to DaVINCI records. See table
    # ORDCOVIDVaccineRecords.
    vax,
    paxlovid,
    # Social
    rurality,
    svi,
    # Clade
    clade,
    # Contact days
    contact_days
  )
missing$n_na <- rowSums(is.na(missing))

# Define exclusions
data_excl <- data |>

  # Join number of missing variables back to each patient/case
  left_join(
    missing |>
      select(patient, index_date, n_na),
    by = c("patient", "index_date")
  ) |>

  # Join flag indicating sequence failed Nextclade quality control
  left_join(
    failed_qc,
    by = join_by(LongAccessionNumberUID == seqName)
  ) |>
  mutate(failed_qc = if_else(is.na(failed_qc), 0L, failed_qc)) |>

  # Compute exclusions
  mutate(
    exclusion = case_when(
      death_flag == 1 ~ "Died within 270 days",
      EverPositive != 1 |
        is.na(EverPositive) |
        is.na(index_date) ~ "No record of COVID positivity",
      PatientCategory != "Veteran" | is.na(PatientCategory) ~
        "Non-veteran or veteran employee",
      legit_sample != 1 | is.na(legit_sample) ~ "No nearby index date",
      n_na > 0 ~ "Missing predictor values",
      outcome == "Digestive" ~ "Gastrointestinal phenotype (too few cases)",
      failed_qc == 1 ~ "Sequence failed Nextclade quality controls",
      TRUE ~ NA_character_
    )
  )
write_feather(data_excl, "data/data_excl.feather")

# Table with info for all samples
samples_all <- data_excl |>
  select(LongAccessionNumberUID, LabSpecimenTakenDateTime, clade_nextstrain)

# Join in indicator for those used in analysis and write to CSV
samples_all2 <- samples_all |>
  left_join(
    data_excl |>
      filter(is.na(exclusion)) |>
      mutate(used = 1L) |>
      select(LongAccessionNumberUID, used),
    by = c("LongAccessionNumberUID")
  ) |>
  mutate(used = if_else(is.na(used), 0L, used))

# Write sample info to CSV
write_excel_csv(
  samples_all2,
  file = "data/sample_data.csv"
)

# Dataset for modeling
data_mod <- data_excl |>
  filter(is.na(exclusion)) |>
  select(
    outcome,
    index_date,
    clade,
    age,
    sex,
    race,
    ethnicity,
    vax,
    rurality,
    svi,
    comorbidity,
    smoking,
    paxlovid,
    dex,
    no_ins,
    contact_days
  ) |>
  mutate(
    # Drop unattested factor levels
    across(where(is.factor), ~ fct_drop(.x)),

    # Log transform contact days, center & scale both contact_days & age
    contact_days = as.numeric(scale(log(contact_days))),
    age = scale(age),

    # Compute week variable from index dates
    week_start = floor_date(index_date, unit = "week"),
    week_index = as.integer(difftime(
      week_start,
      min(week_start, na.rm = TRUE),
      units = "weeks"
    ))
  )

# Extract age scaling variables & write to disk
age_mean <- attr(data_mod$age, "scaled:center")
age_sd <- attr(data_mod$age, "scaled:scale")
write_rds(list(age_mean = age_mean, age_sd = age_sd), "data/age_scale.rds")

# Convert age to numeric to avoid issues with nnet::multinom()
data_mod$age <- as.numeric(data_mod$age)

# Write modeling data to feather
write_feather(data_mod, "data/data_mod.feather")

# Multiple regression models --------------------------------------------------

# Fit with the main predictors we're interested in: viral clade and social variables (age, sex, rurality, svi).
fit1 <- multinom(
  outcome ~ clade + age + sex + rurality + svi,
  data = data_mod,
  MaxNWts = 5000,
  trace = FALSE
)

# Same as fit1, but with additional controls for health (comorbidity, smoking, health care utilization).
fit2 <- multinom(
  outcome ~ clade +
    age +
    sex +
    rurality +
    svi +
    comorbidity +
    smoking +
    contact_days,
  data = data_mod,
  MaxNWts = 5000,
  trace = FALSE
)

# Same as fit2, but with additional controls for treatments (vax & paxlovid).
fit3 <- multinom(
  outcome ~ clade +
    age +
    sex +
    rurality +
    svi +
    comorbidity +
    smoking +
    contact_days +
    vax +
    paxlovid,
  data = data_mod,
  MaxNWts = 5000,
  trace = FALSE
)

# VIF for fit3 predictors. All VIF < 1.75. VIF of 1 is no multicollinearlity; VIF from 1-2 indicates weak multicollinearity.
f_vif <- as.formula(
  y ~ clade +
    age +
    sex +
    rurality +
    svi +
    comorbidity +
    smoking +
    contact_days +
    vax +
    paxlovid
)
vif_lm <- glm(f_vif, data = data_mod |> mutate(y = rnorm(nrow(data_mod))))
vif(vif_lm)

# Risk differences for models 1-3 for target variables clade, age, sex, rurality, svi
rd1 <- avg_comparisons(
  fit1,
  variables = list(
    clade = "pairwise",
    age = 1,
    sex = levels(data_mod$sex),
    rurality = levels(data_mod$rurality),
    svi = levels(data_mod$svi)
  )
) |>
  tidy() |>
  mutate(model = "Base")
rd2 <- avg_comparisons(
  fit2,
  variables = list(
    clade = "pairwise",
    age = 1,
    sex = levels(data_mod$sex),
    rurality = levels(data_mod$rurality),
    svi = levels(data_mod$svi)
  )
) |>
  tidy() |>
  mutate(model = "Base + Health")
rd3 <- avg_comparisons(
  fit3,
  variables = list(
    clade = "pairwise",
    age = 1,
    sex = levels(data_mod$sex),
    rurality = levels(data_mod$rurality),
    svi = levels(data_mod$svi)
  )
) |>
  tidy() |>
  mutate(model = "Base + Health + Treatment")

# Write risk differences to rd.rds
write_rds(
  list(rd1 = rd1, rd2 = rd2, rd3 = rd3),
  "data/rd.rds"
)

# Counterfactual predictions for clade based on fit3
clade_pred <- avg_predictions(
  fit3,
  variables = "clade"
) |>
  tidy() |>
  mutate(
    # Convert estimate & std err to %
    across(c(estimate, std.error), ~ .x * 100),

    # Order clades chronologically
    clade = fct_relevel(clade, "Pre-Delta", "Delta", "Early Omicron"),

    # Order outcomes
    group = fct_relevel(group, "Cardiopulmonary", "Neuropsychiatric")
  )
write_feather(clade_pred, "data/clade_pred.feather")

# Counterfactual predictions for sex based on fit3
sex_pred <- avg_predictions(
  fit3,
  variables = "sex"
) |>
  tidy() |>
  mutate(
    # Convert estimate & std err to %
    across(c(estimate, std.error), ~ .x * 100),

    # Order outcomes
    group = fct_relevel(group, "Cardiopulmonary", "Neuropsychiatric"),

    # Pretty sex labels
    sex = if_else(sex == "F", "Female", "Male")
  )
write_feather(sex_pred, "data/sex_pred.feather")

# Counterfactual predictions for age based on fit3
age_pred <- avg_predictions(
  fit3,
  variables = list(age = seq(min(data_mod$age), max(data_mod$age), by = 0.1))
) |>
  tidy() |>
  mutate(
    # Convert scaled age back to unscaled
    age_unscaled = age * age_sd + age_mean,

    # Convert estimate & std err to %
    across(c(estimate, std.error), ~ .x * 100)
  )
write_feather(age_pred, "data/age_pred.feather")

# Counterfactual predictions for rurality based on fit3
rural_pred <- avg_predictions(
  fit3,
  variables = "rurality"
) |>
  tidy() |>
  mutate(
    # Convert estimate & std err to %
    across(c(estimate, std.error), ~ .x * 100),

    # Order outcomes
    group = fct_relevel(group, "Cardiopulmonary", "Neuropsychiatric")
  )
write_feather(rural_pred, "data/rural_pred.feather")

# Counterfactual predictions for SVI based on fit3
svi_pred <- avg_predictions(
  fit3,
  variables = "svi"
) |>
  tidy() |>
  mutate(
    # Convert estimate & std err to %
    across(c(estimate, std.error), ~ .x * 100),

    # Order outcomes
    group = fct_relevel(group, "Cardiopulmonary", "Neuropsychiatric")
  )
write_feather(svi_pred, "data/svi_pred.feather")

# Unadjusted estimates --------------------------------------------------------

# Target variables
targets <- c("clade", "age", "sex", "rurality", "svi")

# Fit unadjusted model for each target
unadj_rd <- map(targets, function(target) {
  # Create model formula and unadjusted model
  f_unadj <- as.formula(paste0("outcome ~ ", target))
  m_unadj <- multinom(f_unadj, data = data_mod, MaxNWts = 5000, trace = FALSE)

  # Collect risk differences
  if (target == "clade") {
    var_list <- list("pairwise")
    names(var_list) <- target
    comp_unadj <- avg_comparisons(
      m_unadj,
      variables = var_list
    ) |>
      tidy()
  } else {
    comp_unadj <- avg_comparisons(
      m_unadj,
      variables = target
    ) |>
      tidy()
  }
  comp_unadj
}) |>
  list_rbind() |>

  # Clean up table
  mutate(
    term = str_to_title(term),
    term = if_else(term == "Clade", "Variant", term)
  ) |>
  select(
    Term = term,
    Contrast = contrast,
    Outcome = group,
    Estimate = estimate,
    conf.low,
    conf.high
  ) |>
  mutate(across(where(is.numeric), ~ . * 100))
write_feather(unadj_rd, "data/unadj_rd.feather")
