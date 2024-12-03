##---- Kaplan meier curves PACC / no-PACC
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(gridExtra)
library(tidycmprsk)
library(gtsummary)

options(scipen = 0)

survival_km <- function(patient_df, grade=0, label) {
    
    if (!grade %in% c(1,2,3)) {
        grade <- c(1,2,3, NA)
    } else {
        label <- paste0(label,": grade ", grade)   
    }
    
    censor_yr <- 15
    
    surv_df <- patient_df %>% 
        mutate(hist_gra=as.integer(hist_gra),
               PACC_YN = if_else(PACC != 2, 0, 2),
               status_os = if_else(bcd_cmp_ == 1, 1, 0)  # Only death by breast cancer
               ) %>% 
        filter(tumor == 1, hist_gra %in% grade) %>% 
        droplevels()
    
    # Censor the data in 10 years
    surv_df_15 <- surv_df %>%
        mutate(status_os = if_else(surv_dec > censor_yr, 0, status_os),
               surv_dec = if_else(surv_dec > censor_yr, censor_yr, surv_dec))
    
    # Fit model
    fit <- survfit2(
        Surv(surv_dec, status_os) ~ PACC_YN, 
        data = surv_df_15)
    summary(fit)
    
    # Statistics Rho=0 for log-rank, Rho=1 for Peto&Peto mod.
    logrank_test <- survdiff(Surv(surv_dec, status_os) ~ PACC_YN, 
                             data = surv_df_15, rho=0)

    # Prettier plot (can add pvalue of Peto&Peto modf.)
    km_plot <- fit %>%
        ggsurvfit() +
        labs(title = label,
             x = 'Time (years)') +
        add_censor_mark() +
        add_risktable() +
        scale_color_manual(values = c('#619CFF','#F8766D'),
                           labels = c('ECC 0', 'ECC 2')) +
        scale_fill_manual(values = c('#619CFF','#F8766D'),
                          labels = c('ECC 0', 'ECC 2'))+
        ylim(c(0,1)) +
        annotate("text", x = 2.5, y = 0.65,
                 label = paste0("Log-rank p = ", round(logrank_test$pvalue, 4)))
    return(km_plot)
}

recurrence_km <- function(patient_df, grade=0, 
                          rec_type= c("overall", "local","distant"),
                          label) {
    
    rec_type <- match.arg(rec_type)
    
    if (!grade %in% c(1,2,3)) {
            grade <- c(1,2,3, NA)
            label <- paste0(label,": ",rec_type)
    } else {
            label <- paste0(label,": ",rec_type, " - grade", grade)   
        }

    censor_yr <- 10

    rec_df <- patient_df %>%
        filter(tumor == 1, ev1typ20 != 5) %>%
        mutate(hist_gra=as.integer(hist_gra),
               or_rec = if_else(ev1typ20 != 0, 1, 0),
               local_rec = if_else(ev1typ20 == 0, 0, # "rec free"
                                   if_else(ev1typ20 != 4, 1, # "local recurrence"
                                           0)), # "distant recurrence censored"
               dist_rec = if_else(ev1typ20 == 4, 1, # Distant recurrence
                                  0), # rec free and local recurrence censored
               PACC_YN = if_else(PACC != 2, 0, 2),
        ) %>%
        filter(hist_gra %in% grade) %>% 
        droplevels()

    # Censor all the data in 15 years instead of full time.
    rec_df_10 <- rec_df %>%
        mutate(or_rec = if_else(rfs_ny > censor_yr, 0, or_rec),
               local_rec = if_else(rfs_ny > censor_yr, 0, local_rec),
               dist_rec = if_else(rfs_ny > censor_yr, 0, dist_rec),
               rfs_upd = if_else(rfs_ny > censor_yr, censor_yr, rfs_ny))

    if (rec_type == "overall") {
            rec_col <- "or_rec"
    } else if (rec_type == "local") {
            rec_col <- "local_rec"    
    } else {
            rec_col <- "dist_rec"
        } 
    fit_rec <- survfit2(
        Surv(rfs_upd, rec_df_10[[rec_col]]) ~ PACC_YN,
        data = rec_df_10)
    logrank_rec <- survdiff(Surv(rfs_upd, rec_df_10[[rec_col]]) ~ PACC_YN,
                            data = rec_df_10, rho=0)

    km_plot <- fit_rec %>%
        ggsurvfit() +
        labs(title = label,
             x = 'Time (years)',
             y = 'Recurrence-free Remission Probability') +
        add_censor_mark() +
        add_risktable() +
        scale_color_manual(values = c('#619CFF','#F8766D'),
                          labels = c('ECC 0', 'ECC 2')) +
        scale_fill_manual(values = c('#619CFF','#F8766D'),
                         labels = c('ECC 0', 'ECC 2'))+
        ylim(c(0,1)) +
        annotate("text", x = 2.5, y = 0.65,
                 label = paste0("Log-rank p = ", round(logrank_rec$pvalue, 4)))
    return(km_plot)
}

dir <- "~/PACC/TMA/Results/"
tables <- paste0(dir, "ManualScoreTables/")
plot_dir_surv <- paste0(dir, "SurvPlots/OnlyGrade1and2/")
plot_dir_rec <- paste0(dir, "SurvPlots/OnlyGrade1and2/")
patient_data <- read_tsv(paste0(dir,
    "all_patient_unique_with_pacc_full_score_new.csv"), na = "#NULL!") %>%
    mutate(PACC=as.integer(PACC)) %>% 
    filter(exkl == 0, PACC %in% c(0,1,2), tumor==1) %>% 
    droplevels()


# Year cut off to censor
censor_yr <- 15

## Data exploration plots ----
# Distribution of survival
plot1 <- patient_data %>% 
    ggplot(aes(x=surv_dec)) +
    geom_histogram(aes(fill=as.factor(brcadod_)), col='black', alpha=0.5,
                   position = 'stack') +
    labs(title="Distribution of survival time", x="Survival Time (years)", 
         y='Count', fill="") +
    scale_fill_manual(values= c('#7CAE00','#F8766D'), 
                      labels = c('Death by other means', 
                                 "Death by breast cancer")) +
    theme_classic() +
    theme(text = element_text(size=20)); plot1

ggsave(paste0(dir, "ExpPlots/death_by_brca.png"), width = 10, 
       height = 8, 
       units = "in")
## Survival plots ----
for (i in c(0,1,2,3)) {
    surv_plot <- survival_km(patient_df = patient_data, 
                             label="Kaplan-Meier Survival Curve", grade=i)
    filename <- paste0(
        "~/PACC/TMA/Manuscript/KMPlots/survival_km_plot_grade", 
        i, ".tif")        
    
    ggsave(filename, surv_plot, width = 6, height =6, units = "in")
}

## Recurrence plots ----
for (i in c(0,1,2,3)) {
    for (type in c("overall", "local", "distant")) {
            overall <- recurrence_km(patient_df = patient_data, 
                                     rec_type = type,
                                     label="Kaplan-Meier RFR Curve", grade=i)
            filename <- paste0(
                "~/PACC/TMA/Manuscript/KMPlots/", type,
                "_rec_km_plot_grade", 
                i, ".tif")        
            
        ggsave(filename, overall, width = 6, height =6, units = "in")
    }
}
