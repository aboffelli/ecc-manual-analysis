source("renv/activate.R")

suppressPackageStartupMessages(library(tidyverse))
library(survival)
library(broom)      # Tidying statistical results
library(tidycmprsk) # Competing risks

# SCRIPT OVERVIEW
#
# prepare_patient_data: Removes NA values independent variables and sets
# factors, creates a "time" and several "event" columns, the details which vary
# if event_type is survival or recurrence.
#
# fit_all_models: Fits univariable and multivariable coxph (or competing risk)
# models. Fits on entire data set as well as stratifying by grade
#
# Script output is a single large dataframe that contains everything your heart
# could desire. Each line is a 'level' (for numerics, there are no levels in a
# term - so it's just "NA"). Hazard ratios, p-values, confidence intervals (95%)
# as well as 'fit providence' (is it from univariable, multivariable, grade
# stratified, OS, DSS, competing risks, survival, or recurrence). Additionally
# includes the mean and standard deviation of the numeric term, or percent and
# count of the factor term.


# Creates survival/recurrence cols (status_(cr/os/dss)) and censors past the
# cutoff date. Column `time` holds time to event
#
# 0: No event
# 1: Disease specific event
# 2: Non-disease specific event
#
# For OS and DSS, only 0 and 1 exist.
#
# The `local` argument is only used for `event_type = "recurrence"`. It changes
# what is considered to be 'disease' in a DSS sense. If TRUE, local recurrences
# will be considered 'disease', otherwise distant recurrences will.
#
# It probably makes more sense to create many columns prepended with "surv_" and
# "rec_" so this function would only be called once, but here we are.
prepare_patient_data <- function(path,
                                 event_type = c("survival", "recurrence"),
                                 n_years_cutoff = NULL,
                                 local = TRUE) {
  event_type <- match.arg(event_type, c("survival", "recurrence"))

  if (is.null(n_years_cutoff)) {
    if (event_type == "survival") n_years_cutoff <- 15
    if (event_type == "recurrence") n_years_cutoff <- 10
  }

  df <- read_tsv(
    path, na = "#NULL!", col_types = cols(PACC = "i", tumsize = "n")
  ) |>
    filter(
      exkl == 0,
      tumor == 1,
      !is.na(hist_gra),
      !is.na(recgrp),
      !is.na(tumsize)
    ) |>
    mutate(PACC_YN = if_else(PACC == 2, 1, 0)) |>
    filter(!is.na(PACC_YN)) |>
    mutate(
      # read_tsv handles NAs differently than base R, which affects results.
      # Using as.factor for compat.
      across(c(hist_gra, adjbeh, recgrp, rt_treat, PACC_YN, PACC), as.factor)
    )

  # Event-type non-specific manipulations:
  df <- df |>
    # Previously, dropping NAs was done implicitly when filtering for grades to
    # use for analysis.
    #
    # Yes, "NA" is a character string rather than actual NA
    filter(hist_gra != "NA") |>
    droplevels()

  if (event_type == "survival") {
    df <- prepare_survival_df(df, n_years_cutoff)
  }

  if (event_type == "recurrence") {
    df <- prepare_recurrence_df(df, n_years_cutoff, local)
  }

  # tidycmprsk::crr requires event col to be factor, not numeric.
  #
  # Coercing to factor now rather than upstream because factors are slippery and
  # I'd rather deal with them as little as possible.
  df <- mutate(df, status_cr = factor(status_cr))

  df
}
prepare_survival_df <- function(df, n_years_cutoff) {
  df |>
    mutate(
      bcd_cmp_ = if_else(surv_dec > n_years_cutoff, 0, bcd_cmp_),
      time = if_else(surv_dec > n_years_cutoff, n_years_cutoff, surv_dec),
      # Competing risks
      # 0: alive
      # 1: died of breast cancer
      # 2: died of other cause
      status_cr = if_else(bcd_cmp_ %in% 0:1, bcd_cmp_, 2),
      status_os = status_cr != 0,
      status_dss = status_cr == 1
    )
}
prepare_recurrence_df <- function(df, n_years_cutoff, local) {
  local_recurrence_code <- ifelse(local, 1, 2)
  distant_recurrence_code <- ifelse(local, 2, 1)
  df |>
    # No idea what this does, just doing what the old code did
    filter(ev1typ20 != 5) |>
    mutate(
      ev1typ20 = if_else(rfs_ny > n_years_cutoff, 0, ev1typ20),
      time = if_else(rfs_ny > n_years_cutoff, n_years_cutoff, rfs_ny),
      status_cr = case_when(
        ev1typ20 == 0 ~ 0,
        ev1typ20 == 4 ~ distant_recurrence_code,
        .default = local_recurrence_code
      ),
      status_os = status_cr != 0,
      status_dss = status_cr == 1
    )
}

# Fits all models (univariable, multivariable, with overall, disease, and
# competing risk survival). In addition to fitting the entire dataset, also
# stratifies by grade and does the same for each.
#
# As such, model_type (2) * event_type (3) * grades (4) = 24 models will be fit.
# (technically many more since the univariate fits each term seperately)
fit_all_models <- function(df) {
  arg_table <- expand_grid(
    model_type = c("univariable", "multivariable"),
    event_col = c("status_os", "status_dss", "status_cr")
  )

  # Fit without grade stratification
  together <- mapply(
    \(model_type, event_col) fit_coxph(df, model_type, event_col),
    arg_table$model_type,
    arg_table$event_col,
    SIMPLIFY = FALSE
  ) |>
    bind_rows() |>
    mutate(grade_stratum = "all")

  # Fit stratifying by grade
  by_grade <- mapply(
    \(model_type, event_col) {
      tapply(
        df, ~ hist_gra, fit_coxph,
        model_type = model_type,
        event_col = event_col,
        independent_variables = c("age", "rt_treat", "PACC_YN", "recgrp", "tumsize")
      )
    },
    arg_table$model_type, arg_table$event_col, SIMPLIFY = FALSE
  ) |>
    unlist(recursive = FALSE) |>
    bind_rows(.id = "grade_stratum") |>
    # extract grade from 'multivariable.3'
    mutate(grade_stratum = str_extract(grade_stratum, "[^\\.]$"))

  bind_rows(by_grade, together)
}

# Fits a coxph model. Can be univariable or multivariable. Setting `event_col`
# to be `status_cr` with automatically fit with a competing risk model.
fit_coxph <- function(df,
                      model_type = c("univariable", "multivariable"),
                      event_col = "status_dss",
                      independent_variables = c("age", "rt_treat", "PACC_YN", "hist_gra", "recgrp", "tumsize"),
                      time_col = "time") {
  model_type <- match.arg(model_type, c("univariable", "multivariable"))

  if (model_type == "univariable") {
    fit_df <- fit_univariable_coxph(df, time_col, event_col, independent_variables)
  }

  if (model_type == "multivariable") {
    fit_df <- fit_multivariable_coxph(df, time_col, event_col, independent_variables)
  }

  iv_summary <- lapply(independent_variables, \(x) make_iv_summary(df = df, colname = x)) |>
    Reduce(rbind, x = _)

  fit_with_summary <- full_join(fit_df, iv_summary, by = "term")

  # Add some provenance
  fit_with_summary$model_type <- model_type
  fit_with_summary$event_type <- event_col

  fit_with_summary
}

fit_univariable_coxph <- function(df, time_col, event_col, independent_variables) {
  # Use competing risk model if competing risk status is supplied
  cox_fn <- ifelse(event_col == "status_cr", crr, coxph)

  # Make formulas from strings, then fit univariate coxph.
  independent_variables |>
    lapply(
      \(x) {
        paste0("Surv(", time_col, ", ", event_col, ") ~ ", x) |>
          as.formula() |>
          cox_fn(df) |>
          tidy(exponentiate = TRUE, conf.int = TRUE)
      }
    ) |>
    Reduce(rbind, x = _)
}

fit_multivariable_coxph <- function(df, time_col, event_col, independent_variables) {
  formula <- paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste0(independent_variables, collapse = " + ")
  ) |>
    as.formula()

  # Independent variable summary table
  lapply(independent_variables, \(x) {

  })


  cox_fn <- ifelse(event_col == "status_cr", crr, coxph)

  cox_fn(formula, df) |>
    tidy(exponentiate = TRUE, conf.int = TRUE)
}

# Takes in the full dataframe and makes a summary of a given colname.
#
# Produces a dataframe with four columns - term, level, metric_type, metric.
# Metrics are mean and sd if values are numeric, count and percent if factor.
# For numeric metrics, 'level' is NA.
make_iv_summary <- function(df, colname) {
  values <- df[[colname]]
  if (class(values) == "numeric") {
    values_mean <- mean(values, na.rm = TRUE)
    values_sd <- sd(values, na.rm = TRUE)
    value_summary <- data.frame(
      term = colname,
      level = NA, # Dummy column to join with categorical fits
      metric_type = c("mean", "sd"),
      metric = c(values_mean, values_sd)
    )
  }

  if (class(values) == "factor") {
    value_summary <- values |>
      table(useNA = "always") |>
      as.data.frame() |>
      setNames(c("level", "count")) |>
      mutate(
        # term and level get conjoined when fitting factors. This provides a
        # join column to the original fit data.
        term = paste0(colname, level),
        percent = 100 * count / sum(count)
      ) |>
      pivot_longer(
        cols = c("count", "percent"),
        names_to = "metric_type",
        values_to = "metric"
      )
  }
  value_summary
}




recurrence_local <- prepare_patient_data(snakemake@input$patient_data, event_type = "recurrence")
recurrence_distant <- prepare_patient_data(snakemake@input$patient_data, event_type = "recurrence", local = FALSE)
survival <- prepare_patient_data(snakemake@input$patient_data, event_type = "survival")

all_r_local <- fit_all_models(recurrence_local)
all_r_local$dataset <- "recurrence_local"
all_r_distant <- fit_all_models(recurrence_distant)
all_r_distant$dataset <- "recurrence_distant"
all_s <- fit_all_models(survival)
all_s$dataset <- "survival"

all <- bind_rows(all_r_local, all_r_distant, all_s)

write_csv(all, snakemake@output$hr_data)
