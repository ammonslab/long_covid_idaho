source("renv/activate.R")

# Packages ----------------------------------------------------

library(tidyverse)
library(odbc) # Connect to db
library(DBI) # Work with db
library(marginaleffects) # Summarize model results
library(lubridate)
library(gt)
library(DT)
library(car) # vif()
library(consort)
library(feather)
library(ggsignif) # Significance annotations
library(nnet) # Multinomial models
library(janitor) # adorn_totals() to add totals rows
library(furrr) # For future_map()
library(caret) # For confusionMatrix()
library(brglm2) # Penalized multinomial models
library(glmnet) # For glmnet() with family = "multinomial"
library(tictoc) # For timing
library(ggpubr) # Multiple ggplots in same figure
library(ggokabeito)
library(lmtest)
library(gtsummary)
library(flextable)
library(english)
library(tidygraph)
library(ggraph)
library(patchwork) # Needed to align plots in a panel figure
library(here) # Knit from directories other than project
library(modelsummary) # Multi-model tables
library(splines)

# Resolve conflicts
select <- dplyr::select
filter <- dplyr::filter

# ggplot theme -----------------------------------------

theme_set(theme_minimal())

# R-squared ---------------------------------------------------

rsq <- function(y_actual, y_predicted) {
  cor(y_actual, y_predicted)^2
}

# Logistic regression models & summary ------------------------

fit_sum <- function(data, outcome, predictors) {
  # User message
  message(paste0("Modeling ", outcome, "..."))

  # Vector of possible outcomes
  outcomes <- c("long_covid_flag", "cardio_flag", "neuro_flag", "dig_flag")

  # Outcome cols to drop
  drop <- outcomes[-which(outcomes == outcome)]

  # Modeling dataset
  data <- data %>%
    select(all_of(c(outcome, predictors)))

  # Model formula
  formula <- as.formula(paste(outcome, "== 1 ~ ."))

  # Fit model
  fit <- glm(formula, data = data, family = binomial(link = "logit"))

  # Compute r-squared
  rsq <- rsq(data %>% pull(outcome), fitted(fit))

  # Vector of predictors that don't require pairwise comparisons
  non_pair_preds <- names(data)[
    -which(
      names(data) %in%
        c(outcome, "clade_nextstrain", "Race", "Ethnicity", "cdc")
    )
  ]

  # Population-averaged adjusted risk differences
  diff <- bind_rows(
    avg_comparisons(fit, variables = non_pair_preds) %>%
      tidy(),
    avg_comparisons(fit, variables = list(Race = "pairwise")) %>%
      tidy(),
    avg_comparisons(fit, variables = list(Ethnicity = "pairwise")) %>%
      tidy(),
    avg_comparisons(fit, variables = list(cdc = "pairwise")) %>%
      tidy(),
    avg_comparisons(fit, variables = list(clade_nextstrain = "pairwise")) %>%
      tidy()
  )

  # Apply Bonferroni correction
  diff2 <- diff %>%
    mutate(statsig = if_else(p.value < 0.05 / nrow(diff), 1L, 0L))

  # Combine into list
  results <- list(fit = fit, rsq = rsq, diff = diff2)
}

# Summarize columns in a dataframe ------------------------

sum_df <- function(data) {
  sum <- purrr::map(names(data), function(var) {
    n_missing <- sum(is.na(data %>% pull({{ var }})))

    if (is.numeric(data %>% pull({{ var }}))) {
      min <- min(data %>% pull({{ var }}), na.rm = TRUE)
      max <- max(data %>% pull({{ var }}), na.rm = TRUE)
      mean <- mean(data %>% pull({{ var }}), na.rm = TRUE)
      n_zeroes <- sum(data %>% pull({{ var }}) == 0, na.rm = TRUE)
    } else {
      min <- NA
      mean <- NA
      max <- NA
      n_zeroes <- NA
    }
    tibble(
      variable = {{ var }},
      n_missing = n_missing,
      min = min,
      max = max,
      mean = mean,
      n_zeroes = n_zeroes
    ) %>%
      mutate(
        p_missing = n_missing / nrow(data),
        p_zeroes = n_zeroes / nrow(data)
      ) %>%
      select(variable, p_missing, n_missing:p_zeroes)
  }) %>%
    list_rbind() %>%
    arrange(desc(p_missing))
  return(sum)
}

# Test for differences in outcomes by clade ------------------------

# Function that takes a model and clade as input and computes within-clade
# adjusted population-averaged outcome contrasts.
clade_contrasts <- function(model, clade) {
  message(paste0("Processing ", clade, "..."))
  avg_predictions(
    model,
    type = "probs",
    variables = list(clade_nextstrain = clade),
    hypothesis = "pairwise"
  ) |>
    tidy() |>
    mutate(clade = clade) |>
    select(clade, term:conf.high) %>%
    return()
}

# Plot predictions for all phenotypes -----------------------------------------

# Function that takes predictions and a variable as input & plots all phenotype
# predictions for the variable in a single panel.
plot_pheno_all <- function(pred, var) {
  # Figure parameters
  pd <- position_dodge(0.3)
  error_width <- 0.4
  error_linewidth <- 0.2
  point_size <- 2.1
  point_alpha <- 0.5

  # Different mutations on data depending on variable
  if (var == "cdc") {
    fig_data <- pred |>
      filter(term == var, group != "no_lc") |>
      mutate(
        x = str_wrap(x, 6),
        x = if_else(x == "Non-Core", "Non-\nCore", x),
        x = fct_relevel(
          x,
          "Large\nMetro",
          "Medium\nMetro",
          "Small\nMetro",
          "Micropolitan",
          "Non-\nCore"
        )
      )
  } else if (var == "svi") {
    fig_data <- pred |>
      filter(term == var, group != "no_lc") |>
      mutate(x = str_to_title(x), x = fct_relevel(x, "Low", "High"))
  } else if (var == "clade_nextstrain") {
    fig_data <- pred |>
      filter(term == var, group != "no_lc") |>
      mutate(
        x = fct_relevel(
          x,
          "20A lineage",
          "20B lineage",
          "Delta",
          "Omicron",
          "Omicron descendant",
          "Recombinant"
        )
      )
  } else {
    fig_data <- pred |>
      filter(term == var, group != "no_lc") |>
      mutate(x = str_to_title(x), x = str_wrap(x, 22))
  }
  fig <- fig_data |>
    ggplot(aes(
      x = x,
      y = estimate,
      color = group,
      group = group
    )) +
    geom_errorbar(
      aes(ymin = estimate - std.error, ymax = estimate + std.error),
      position = pd,
      linewidth = error_linewidth,
      width = error_width
    ) +
    geom_line(position = pd) +
    geom_point(position = pd, size = point_size, alpha = point_alpha) +
    scale_color_okabe_ito(
      breaks = c("Neuropsychiatric", "Cardiopulmonary", "Digestive"),
      labels = c("Neuropsychiatric", "Cardiopulmonary", "Digestive")
    ) +
    labs(color = "Phenotype", x = var, y = "") +
    theme(
      panel.grid.major = element_line(color = "gray80", linewidth = 0.1),
      legend.position = "top"
    )
  return(fig)
}

# Plot predictions for individual phenotype -----------------------------------

# Function that takes a table of statsig comparisons, a table of predictions for
# a predictor, and a table specifying the order of x-axis labels and outputs
# prediction figures for each phenotype, with annotations indicating statsig
# comparisons.
plot_pheno_indiv <- function(comp, pred, order) {
  # Figure parameters
  pd <- position_dodge(0.3)
  error_width <- 0.4
  error_linewidth <- 0.2
  point_size <- 2.1
  point_alpha <- 0.5

  # Extract term
  term <- comp |>
    pull(term) |>
    unique()

  # Statsig annotations for clade
  anno <- comp |>
    select(group, contrast) |>
    separate_wider_delim(cols = contrast, delim = " - ", names = c("a", "b")) |>
    mutate(across(a:b, ~ fct_relevel(.x, order$label))) |>
    left_join(order |> rename(a_order = order), by = join_by(a == label)) |>
    left_join(order |> rename(b_order = order), by = join_by(b == label)) |>
    mutate(
      start = if_else(a_order < b_order, a, b),
      end = if_else(a_order > b_order, a, b)
    ) |>
    select(group, start, end)

  # Make figure for each phenotype
  pheno_figs <- map(
    c("Neuropsychiatric", "Cardiopulmonary", "Digestive"),
    function(pheno) {
      # Define color to be used for phenotype
      pheno_color <- case_when(
        pheno == "Neuropsychiatric" ~ 3,
        pheno == "Cardiopulmonary" ~ 1,
        pheno == "Digestive" ~ 2,
        TRUE ~ 999
      )

      # Get data for phenotype
      pheno_data <- pred |>
        filter(group == pheno) |>
        mutate(upper = estimate + std.error)

      # Upper defines the upper bound where data is plotted. We'll put statsig
      # annotations above this.
      upper <- max(pheno_data$upper) + 1

      # Annotations to add for phenotype
      pheno_anno <- anno |>
        filter(group == pheno)

      # Create base figure
      base_fig <- pheno_data |>
        ggplot(aes(
          x = x,
          y = estimate,
          color = group,
          group = group
        )) +
        geom_errorbar(
          aes(ymin = estimate - std.error, ymax = estimate + std.error),
          position = pd,
          linewidth = error_linewidth,
          width = 0.2
        ) +
        geom_line(position = pd) +
        geom_point(position = pd, size = point_size, alpha = point_alpha) +
        scale_color_okabe_ito(order = pheno_color) +
        labs(
          x = term,
          y = "Phenotype (%)",
          color = "Phenotype"
        ) +
        theme(
          panel.grid.major = element_line(color = "gray80", linewidth = 0.1),
          legend.position = "top"
        )

      # Add annotations to base figure if it has any
      if (nrow(pheno_anno) > 0) {
        # Table with all data needed for annotation
        signif_tbl <- pheno_anno |>
          mutate(y = upper + seq(0, nrow(pheno_anno) - 1, 1), label = "")

        # Add annotations
        base_fig <- base_fig +
          geom_signif(
            data = signif_tbl,
            aes(
              xmin = start,
              xmax = end,
              y_position = y,
              annotations = label,
              group = group
            ),
            manual = TRUE,
            color = "black",
            size = 0.2
          )
      }

      # Return this
      base_fig
    }
  )

  return(pheno_figs)
}

# Fx to do quantile binning ---------------------------------------------------
qbin <- function(x, q = 5) {
  cut(
    x,
    breaks = quantile(x, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE),
    include.lowest = TRUE,
    ordered_result = TRUE
  )
}

# Fx to check for factor x factor support -------------------------------------
support_ff <- function(
  df,
  v1,
  v2,
  outcome = "outcome",
  min_n = 8,
  per_class = FALSE
) {
  if (per_class) {
    tab <- xtabs(~ df[[v1]] + df[[v2]] + df[[outcome]])
    ok <- all(tab >= min_n)
  } else {
    tab <- xtabs(~ df[[v1]] + df[[v2]])
    ok <- all(tab >= min_n)
  }
  list(ok = ok, table = tab)
}
# Fx to check for factor x continuous support ---------------------------------
support_fc <- function(
  df,
  f,
  x,
  outcome = "outcome",
  min_n = 8,
  q = 5,
  per_class = FALSE
) {
  xb <- qbin(df[[x]], q = q)
  if (per_class) {
    tab <- xtabs(~ df[[f]] + xb + df[[outcome]])
    ok <- all(tab >= min_n)
  } else {
    tab <- xtabs(~ df[[f]] + xb)
    ok <- all(tab >= min_n)
  }
  list(ok = ok, table = tab)
}

# Fx to check for continuous x continuous support -----------------------------
support_cc <- function(
  df,
  x1,
  x2,
  outcome = "outcome",
  min_n = 8,
  q = 5,
  per_class = FALSE
) {
  x1b <- qbin(df[[x1]], q = q)
  x2b <- qbin(df[[x2]], q = q)
  if (per_class) {
    tab <- xtabs(~ x1b + x2b + df[[outcome]])
    ok <- all(tab >= min_n)
  } else {
    tab <- xtabs(~ x1b + x2b)
    ok <- all(tab >= min_n)
  }
  list(ok = ok, table = tab)
}

# Fx to check support overall interaction support -----------------------------
keep_supported <- function(
  df,
  keep,
  outcome = "outcome",
  min_n = 8,
  q = 5,
  per_class = FALSE
) {
  is_interaction <- grepl(":", keep)
  inter_terms <- keep[is_interaction]

  ok_terms <- c()
  for (trm in inter_terms) {
    parts <- strsplit(trm, ":")[[1]]

    # Map back to original variable names
    base_of <- function(s) {
      str_extract(s, paste0(predictors, collapse = "|"))
    }
    v1 <- base_of(parts[1])
    v2 <- base_of(parts[2])

    # Decide which support fx to use
    type1 <- if (is.numeric(df[[v1]])) "cont" else "fact"
    type2 <- if (is.numeric(df[[v2]])) "cont" else "fact"

    ok <- switch(
      paste(type1, type2, sep = "_"),
      fact_fact = support_ff(df, v1, v2, outcome, min_n, per_class)$ok,
      fact_cont = support_fc(df, v1, v2, outcome, min_n, q, per_class)$ok,
      cont_fact = support_fc(df, v2, v1, outcome, min_n, q, per_class)$ok,
      cont_cont = support_cc(df, v1, v2, outcome, min_n, q, per_class)$ok
    )
    if (isTRUE(ok)) ok_terms <- c(ok_terms, paste0(v1, ":", v2))
  }
  c(keep[!is_interaction], ok_terms)
}

# Fx to compute pairwise adjusted risk diff -----------------------------------
pairwise_comp <- function(model, term) {
  message(paste("Computing pairwise comparisons for", term))

  # Create list to feed to avg_comparisons() variable argument
  term_list <- list("pairwise")
  names(term_list) <- term

  # Compute adjusted risk differences
  comp <- avg_comparisons(
    model,
    variables = term_list
  ) |>
    tidy()
  return(comp)
}
