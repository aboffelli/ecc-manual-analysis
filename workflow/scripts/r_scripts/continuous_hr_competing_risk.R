source("renv/activate.R")

suppressPackageStartupMessages(library(tidyverse))
library(survival)
library(finalfit)

surv_hr <- function(patient_df, grade = c(1, 2, 3), label, censor_yr = 15) {
  stopifnot(all(grade %in% c(1, 2, 3)))

  explanatory <- c("age", "rt_treat", "PACC_YN", "hist_gra", "recgrp", "tumsize")
  os_label <- "Overall Survival"

  if (length(grade) == 1) {
    # Since one grade is selected, remove grade (which is now constant) from the
    # regression
    explanatory <- explanatory[-which(explanatory == "hist_gra")]
    os_label <- paste0(os_label, ": NHG ", grade)
  }

  surv_df_15 <- patient_df %>%
    mutate(
      bcd_cmp_ = if_else(surv_dec > censor_yr, 0, bcd_cmp_),
      surv_dec = if_else(surv_dec > censor_yr, censor_yr, surv_dec)
    ) %>%
    filter(hist_gra %in% grade) |>
    mutate(across(c(recgrp, hist_gra), droplevels))

  surv_df_15 <- surv_df_15 |>
    mutate(
      # Competing risks regression
      # 0: alive
      # 1: died of breast cancer
      # 2: died of other cause
      status_crr = if_else(bcd_cmp_ %in% 0:1, bcd_cmp_, 2),
      status_os = status_crr != 0,
      status_dss = status_crr == 1,

      # Rename the variables (for multivariate analysis).
      age = ff_label(age, "Age (years)"),
      rt_treat = rt_treat |>
        fct_recode("No" = "0", "Yes" = "1") |>
        ff_label("Radiotherapy"),
      hist_gra = ff_label(hist_gra, "NHG"),
      recgrp = recgrp |>
        fct_recode(
          "Luminal A" = "1",
          "Luminal B" = "2",
          "HER2" = "3",
          "TN" = "4"
        ) |>
        fct_relevel("HER2", "Luminal A", "Luminal B", "TN") |>
        ff_label("Cancer subtype"),
      tumsize = ff_label(tumsize, "Tumor size"),
      PACC_YN = PACC_YN |>
        fct_recode("-" = "0", "+" = "1") |>
        ff_label("ECC Score")
    )
  # Survival objects
  dependent_os <- "Surv(rfs_upd, status_os)"
  dependent_dss <- "Surv(surv_dec, status_dss)"
  dependent_crr <- "Surv(surv_dec, status_crr)"

  os_cph <- surv_df_15 %>%
    finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
    rename(!!os_label := label) %>%
    rename(" " = levels) %>%
    rename("  " = all)

  my_table <- make_table(
    surv_df_15, explanatory, dependent_dss, dependent_crr, label
  )

  return(list(os_cph, my_table))
}

rec_hr <- function(patient_df, grade = c(1, 2, 3), label, local = TRUE, censor_yr = 10) {
  stopifnot(all(grade %in% c(1, 2, 3)))

  explanatory <- c("age", "rt_treat", "PACC_YN", "hist_gra", "recgrp", "tumsize")
  or_label <- "Overall Recurrence"

  if (length(grade) == 1) {
    explanatory <- explanatory[-which(explanatory == "hist_gra")]
    # NOTE: extra space ("NHG  ") added for compatibility with old file
    or_label <- paste0(or_label, ": NHG  ", grade)
  }

  rec_df <- patient_df %>%
    filter(ev1typ20 != 5)

  rec_df_10 <- rec_df %>%
    mutate(ev1typ20 = if_else(rfs_ny > censor_yr, 0, ev1typ20),
           rfs_upd = if_else(rfs_ny > censor_yr, censor_yr, rfs_ny)) %>%
    filter(hist_gra %in% grade) |>
    mutate(across(c(recgrp, hist_gra), droplevels))

  if (local) {
    rec_df_10 <- rec_df_10 |>
      mutate(
        # Competing risks regression
        # 0: recurrence free
        # 1: local recurrence
        # 2: distant recurrence
        status_crr = case_when(
          ev1typ20 == 0 ~ 0,
          ev1typ20 == 4 ~ 2,
          .default = 1
        )
      )
  } else {
    rec_df_10 <- rec_df_10 |>
      mutate(
        # Competing risks regression
        # 0: recurrence free
        # 1: distant recurrence
        # 2: local recurrence
        status_crr = case_when(
          ev1typ20 == 0 ~ 0,
          ev1typ20 == 4 ~ 1,
          .default = 2
        )
      )
  }

  rec_df_10 <- rec_df_10 |>
    mutate(
      status_os = status_crr != 0,
      status_dss = status_crr == 1
    )

  rec_df_10 <- rec_df_10 %>%
    mutate(
      # Rename the variables.
      age = ff_label(age, "Age (years)"),

      rt_treat = rt_treat %>%
        fct_recode("No" = "0",
                   "Yes" = "1") %>%
        ff_label("Radiotherapy"),

      hist_gra = ff_label(hist_gra, "NHG"),

      recgrp = recgrp %>%
        fct_recode(
          "Luminal A" = "1",
          "Luminal B" = "2",
          "HER2" = "3",
          "TN" = "4"
        ) %>%
        ff_label("Cancer subtype"),
      tumsize = ff_label(tumsize, "Tumor size"),
      PACC_YN = PACC_YN %>%
        fct_recode("-" = "0",
                   "+" = "1") %>%
        ff_label("ECC Score")
    )


  # Set the multivariable columns and the tests
  dependent_os <- "Surv(rfs_upd, status_os)"
  dependent_dss <- "Surv(rfs_upd, status_dss)"
  dependent_crr <- "Surv(rfs_upd, status_crr)"

  or_cph <- rec_df_10 %>%
    finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
    rename(!!or_label := label) %>%
    rename(" " = levels) %>%
    rename("  " = all)

  my_table <- make_table(
    rec_df_10, explanatory, dependent_dss, dependent_crr, label
  )

  return(list(or_cph, my_table))
}

make_table <- function(df, explanatory, dependent_dss, dependent_crr, label) {
  # Generate the table by merging all regressions together.
  df |>
    # Summary table
    summary_factorlist(dependent_dss, explanatory,
                       column = TRUE, fit_id = TRUE) |>
    # CPH univariable
    ff_merge(
      explanatory |>
        map(~ df |>
              coxph(as.formula(paste(dependent_dss, .x, sep = " ~ ")),
                    data = _) |>
              fit2df(estimate_suffix = " CPH univariable")) |>
        bind_rows()
    ) |>
    # CPH multivariable
    ff_merge(
      df |>
        coxphmulti(dependent_dss, explanatory) |>
        fit2df(estimate_suffix = " CPH multivariable")
    ) |>
    # Fine and Gray competing risks regression
    ff_merge(
      df |>
        crrmulti(dependent_crr, explanatory) |>
        fit2df(estimate_suffix = " competing risks multivariable")
    ) |>
    select(-fit_id, -index) |>
    dependent_label(df, label)
}

if (!interactive()) {
  input_file <- snakemake@input$patient_data
  o <- snakemake@output
  out_surv <- o$hr_surv
  out_surv_grades <- c(o$hr_surv_grade_1, o$hr_surv_grade_2, o$hr_surv_grade_3)
  out_rec <- o$hr_rec
  out_rec_grades <- c(o$hr_rec_grade_1, o$hr_rec_grade_2, o$hr_rec_grade_3)
}

if (interactive()) {
  input_file <- "./results/all_patient_unique_with_ecc_hif1a_gex.txt"
  out_surv <- "./results/hr/cont_os_with_grade.csv"
  out_surv_grades <- paste0("./results/hr/cont_os_with_grade", 1:3, ".csv")
  out_rec <- "./results/hr/cont_or_with_grade.csv"
  out_rec_grades <- paste0("./results/hr/cont_or_with_grade", 1:3, ".csv")
}

options(scipen = 0)

patient_data <- read_tsv(
  input_file, na = "#NULL!", col_types = cols(PACC = "i", tumsize = "n")
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
    # read_tsv handles NAs differently than base R, which affects results. Using
    # as.factor for compat.
    across(c(hist_gra, adjbeh, recgrp, rt_treat, PACC_YN, PACC), as.factor)
  )

## Survival ----
# Survival with all grades
all_grades <- surv_hr(patient_df = patient_data,
                      label = "Survival")

write_tsv(all_grades[[2]], out_surv)

# Survival for all grades separately
for (i in c(1, 2, 3)) {
  new_table <- surv_hr(patient_df = patient_data,
                       grade = i,
                       label = paste("Survival: Grade", i))
  write_tsv(new_table[[2]], out_surv_grades[i])
}

## Overall Recurrence ----
# Recurrence with all grades
all_grades <- rec_hr(patient_data,
                     label = "Local recurrence",
                     local = TRUE)

write_tsv(all_grades[[1]], out_rec)

# Recurrence for all grades separately
for (i in c(1, 2, 3)) {
  new_table <- rec_hr(patient_data, grade = i,
                      label = paste("Local recurrence: Grade", i),
                      local = TRUE)

  write_tsv(new_table[[1]], out_rec_grades[i])
}
