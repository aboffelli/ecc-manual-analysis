#---- Kaplan meier curves PACC / no-PACC
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(gridExtra)
library(tidycmprsk)
library(gtsummary)

options(scipen = 0)

dir <- "~/PACC/TMA/Results/"
tables <- paste0(dir, "ManualScoreTables/")
plot_dir_surv <- paste0(dir, "SurvPlots/OnlyGrade1and2/")
plot_dir_rec <- paste0(dir, "SurvPlots/OnlyGrade1and2/")
patient_data <- read_tsv(
    paste0(dir,
           "all_patients_with_TMA_PACC_HIF.csv")) %>%
    mutate(PACC=as.integer(PACC)) %>% 
    filter(PACC %in% c(0,1,2)) %>% 
    droplevels() %>% 
    mutate(PACC_YN = if_else(PACC==0, 0, 1),
           age_gap = if_else(age >= 50, ">=50", "<50"),
           tumor_size = if_else(tumsize <= 20, "<=20mm", ">20mm"))

    
# All ----
stratified_summary <- patient_data %>%
    group_by(PACC) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

## Age ----
group_summary <- patient_data %>%
    filter(!is.na(age_gap)) %>%
    group_by(age_gap) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(age_gap)) %>%
    group_by(age_gap, PACC) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$age, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Tumour size ----
group_summary <- patient_data %>%
    filter(!is.na(tumor_size)) %>% 
    group_by(tumor_size) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(tumor_size)) %>% 
    group_by(tumor_size, PACC) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$tumsize, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Ki67 ----
group_summary <- patient_data %>%
    filter(!is.na(ki67_pos)) %>% 
    group_by(ki67_pos) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(ki67_pos)) %>% 
    group_by(ki67_pos, PACC) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$ki67_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=spearman_test$p.value)
print(rho)
print(p_value)
## NHG ----
group_summary <- patient_data %>%
    filter(!is.na(hist_gra)) %>% 
    group_by(hist_gra) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(hist_gra)) %>% 
    group_by(hist_gra, PACC) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$hist_gra, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)
## ER ----
group_summary <- patient_data %>%
    filter(!is.na(er_pos)) %>% 
    group_by(er_pos) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(er_pos)) %>% 
    group_by(er_pos, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$er_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)
## PR ----
group_summary <- patient_data %>%
    filter(!is.na(pgr_pos)) %>% 
    group_by(pgr_pos) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(pgr_pos)) %>% 
    group_by(pgr_pos, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$pgr_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## HER2 ----
group_summary <- patient_data %>%
    filter(!is.na(her2_com)) %>% 
    group_by(her2_com) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(her2_com)) %>% 
    group_by(her2_com, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$her2_com, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## HIF ----
group_summary <- patient_data %>%
    filter(!is.na(HIF1A_ACTIVATION)) %>% 
    group_by(HIF1A_ACTIVATION) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(HIF1A_ACTIVATION)) %>% 
    group_by(HIF1A_ACTIVATION, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$HIF1A_ACTIVATION, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Grade markers - atypia ----
group_summary <- patient_data %>%
    filter(!is.na(atypia)) %>% 
    group_by(atypia) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(atypia)) %>% 
    group_by(atypia, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$atypia, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Grade markers - mitosis ----
group_summary <- patient_data %>%
    filter(!is.na(mitosis)) %>% 
    group_by(mitosis) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(mitosis)) %>% 
    group_by(mitosis, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$mitosis, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Grade markers - tubular formation ----
group_summary <- patient_data %>%
    filter(!is.na(tubuli)) %>% 
    group_by(tubuli) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(tubuli)) %>% 
    group_by(tubuli, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$tubuli, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Adjuvant therapy ----
group_summary <- patient_data %>%
    filter(!is.na(adjbeh)) %>% 
    group_by(adjbeh) %>%
    summarise(
        n = n(),
        percentage = n() / nrow(patient_data) * 100
    )
print(group_summary)

stratified_summary <- patient_data %>%
    filter(!is.na(adjbeh)) %>% 
    group_by(adjbeh, PACC) %>%
    summarise(
        n = n(),
        percentage = n / nrow(patient_data) * 100
    ) %>%
    ungroup()
print(stratified_summary)

spearman_test <- cor.test(patient_data$adjbeh, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)
