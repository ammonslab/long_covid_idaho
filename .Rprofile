source("renv/activate.R")

# Packages ----------------------------------------------------

library(pak)
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

# Custom ggplot theme
my_theme <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_line(
        color = "gray50",
        linewidth = 0.01
      ),
      strip.background = element_blank()
    )
}
theme_set(my_theme())
