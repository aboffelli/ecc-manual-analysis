##---- Competing risk regression
library(tidyverse)
library(cmprsk)
library(finalfit)
library(gt)

custom_theme <- function(){theme_minimal() +
    theme(
        # Remove y-axis line
        axis.line.x = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14)
        
        
    )}

surv_hr <- function(patient_df, grade=0, label) {
    # Select the variable to be included in the regression
    if (grade %in% c(1,2,3)) {  # If one grade is selected, remove it from the regression
        explanatory   <- c("age", 
                           "rt_treat", 
                           "PACC_YN", 
                           "recgrp", 
                           "tumsize")
        os_label <- paste("Overall Survival: NHG ", grade)
    } else {  # Include grade in the regression if using all grades
        grade = c(1,2,3)
        explanatory   <- c("age", 
                           "rt_treat",
                           "PACC_YN", 
                           "hist_gra", 
                           "recgrp", 
                           "tumsize")
        os_label <- paste("Overall Survival")
    }
    
    censor_yr <- 15
    surv_df <- patient_df %>% 
        filter(PACC_YN %in% c(0,1,2)) %>% 
        mutate(rt_treat=as.factor(rt_treat), 
               PACC_YN=as.factor(PACC_YN),
               PACC=as.factor(PACC)) %>% 
        droplevels()
    
    # Censor the data based on the cut off year (15).
    surv_df_15 <- surv_df %>%
        mutate(
            bcd_cmp_ = if_else(surv_dec > censor_yr, 0, bcd_cmp_),
            surv_dec = if_else(surv_dec > censor_yr, censor_yr, surv_dec),
        ) %>% 
        filter(!is.na(hist_gra),  # Filter these out if using multivariate.
               !is.na(recgrp),
               !is.na(tumsize),
               hist_gra %in% grade) %>%
        droplevels()
    
    # Recode the data 
    # For survival, the competing risk is important.
    surv_df_15 <- surv_df_15 %>% 
        mutate(
            # Overall survival
            status_os = if_else(bcd_cmp_ == 0, 0, # "still alive"
                                1), # "died of bc" or "died of other causes"
            
            # Diease-specific survival
            status_dss = if_else(bcd_cmp_ == 0, 0, # "still alive"
                                 if_else(bcd_cmp_ == 1, 1, # "died of bc"
                                         0)), # "died of other causes is censored"
            
            # Competing risks regression
            status_crr = if_else(bcd_cmp_ == 0, 0, # "still alive"
                                 if_else(bcd_cmp_ == 1, 1, # "died of bc"
                                         2)), # "died of other causes"
            
            # Rename the variables (for multivariate analysis).
            age = ff_label(age, "Age (years)"),
            
            rt_treat= factor(rt_treat) %>% 
                fct_recode("No" = "0",
                           "Yes" = "1") %>% 
                ff_label("Radiotherapy"),
            
            hist_gra=factor(hist_gra) %>%
                ff_label("NHG"),
            
            recgrp=factor(recgrp) %>%
                fct_recode(
                    "Luminal A" = "1",
                    "Luminal B" = "2",
                    "HER2" = "3",
                    "TN" = "4"
                ) %>%
                ff_label("Cancer subtype"),
            tumsize=ff_label(as.numeric(tumsize), "Tumor size"),
            PACC_YN=factor(PACC_YN) %>% 
                fct_recode("-"= "0",
                           "+" = "1") %>% 
                ff_label("ECC Score")
        )
    
    print(table(surv_df_15$PACC_YN, surv_df_15$status_crr))
    print(table(surv_df_15$hist_gra, surv_df_15$status_crr))
    # Survival objects
    dependent_os <- "Surv(rfs_upd, status_os)"
    dependent_dss <- "Surv(surv_dec, status_dss)"
    dependent_crr <- "Surv(surv_dec, status_crr)"
    
    os_cph <- surv_df_15 %>% 
        finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
        rename(!!os_label := label) %>% 
        rename(" " = levels) %>% 
        rename("  " = all)
    
    # Generate the table by merging all regressions together.
    my_table <- surv_df_15 %>%
        # Summary table
        summary_factorlist(dependent_dss, explanatory, 
                           column = TRUE, fit_id = TRUE) %>%
        # CPH univariable
        ff_merge(
            explanatory %>%
                map(~ surv_df_15 %>%
                        coxph(as.formula(paste(dependent_dss, .x, sep = " ~ ")), 
                              data = .) %>%
                        fit2df(estimate_suffix = " CPH univariable")
                )%>%
                bind_rows()
        ) %>%
        # CPH multivariable
        ff_merge(
            surv_df_15 %>%
                coxphmulti(dependent_dss, explanatory) %>%
                fit2df(estimate_suffix = " CPH multivariable")
        ) %>%
        # Fine and Gray competing risks regression
        ff_merge(
            surv_df_15 %>%
                crrmulti(dependent_crr, explanatory) %>%
                fit2df(estimate_suffix = " competing risks multivariable")
        ) %>%
        select(-fit_id,-index) %>%
        dependent_label(surv_df_15, label)
    
    return(list(os_cph,my_table))
    
}

rec_hr <- function(patient_df, grade=0, label, local=TRUE) {
    if (grade %in% c(1,2,3)) {
        explanatory   <- c("age", 
                           "rt_treat", 
                           "PACC_YN", 
                           "recgrp", 
                           "tumsize")
        or_label <- paste("Overall Recurrence: NHG ", grade)
    } else {
        grade = c(1,2,3)
        explanatory   <- c("age",
                           "rt_treat", 
                           "PACC_YN", 
                           "hist_gra", 
                           "recgrp", 
                           "tumsize")
        or_label <- "Overall Recurrence"
    }
    
    censor_yr <- 10  # Rec we use 10 years censor 
    rec_df <- patient_df %>% 
        filter(PACC_YN %in% c(0,1,2), tumor == 1,
               ev1typ20 != 5) %>% 
        mutate(rt_treat=as.factor(rt_treat),
               PACC_YN=as.factor(PACC_YN)) %>% 
        droplevels()
    
    rec_df_10 <- rec_df %>%
        mutate(ev1typ20 = if_else(rfs_ny > censor_yr, 0, ev1typ20),
               rfs_upd = if_else(rfs_ny > censor_yr, censor_yr, rfs_ny)) %>% 
        filter(!is.na(hist_gra),  # Filter these out if using multivariate.
               !is.na(recgrp),
               !is.na(tumsize),
               hist_gra %in% grade
               ) %>%
        droplevels()
   
    rec_df_10 <- rec_df_10 %>% 
        mutate(
            # Rename the variables.
            age = ff_label(age, "Age (years)"),
            
            rt_treat= factor(rt_treat) %>% 
                fct_recode("No" = "0",
                           "Yes" = "1") %>% 
                ff_label("Radiotherapy"),
            
            hist_gra=factor(hist_gra) %>%
                ff_label("NHG"),
            
            recgrp=factor(recgrp) %>%
                fct_recode(
                    "Luminal A" = "1",
                    "Luminal B" = "2",
                    "HER2" = "3",
                    "TN" = "4"
                ) %>%
                ff_label("Cancer subtype"),
            tumsize=ff_label(as.numeric(tumsize), "Tumor size"),
            PACC_YN=factor(PACC_YN) %>% 
                fct_recode("-"= "0",
                           "+" = "1") %>% 
                ff_label("ECC Score")
        )
    
    if (local == TRUE) {
        rec_df_10 <- rec_df_10 %>% 
            mutate(
            # Overall survival
            status_os = if_else(ev1typ20 == 0, 0, # "rec free"
                                1), # "recurrence any type"
            
            # Diease-specific survival
            status_dss = if_else(ev1typ20 == 0, 0, # "rec free"
                                 if_else(ev1typ20 != 4, 1, # "local recurrence"
                                         0)), # "distant recurrence censored"
            
            # Competing risks regression
            status_crr = if_else(ev1typ20 == 0, 0, # "rec free"
                                 if_else(ev1typ20 != 4, 1, # local recurrence"
                                         2)), # "distant recurrence"
        )
        
    } else {
        rec_df_10 <- rec_df_10 %>% 
            mutate(
            # Overall survival
            status_os = if_else(ev1typ20 == 0, 0, # "rec free"
                                1), # "recurrence any type"
            
            # Diease-specific survival
            status_dss = if_else(ev1typ20 == 0, 0, # "rec free"
                                 if_else(ev1typ20 == 4, 1, # "distant recurrence"
                                         0)), # "local recurrence censored"
            
            # Competing risks regression
            status_crr = if_else(ev1typ20 == 0, 0, # "rec free"
                                 if_else(ev1typ20 == 4, 1, # distant recurrence"
                                         2)), # "local recurrences"
        )
    }
    
    print(table(rec_df_10$PACC_YN, rec_df_10$status_crr))
    print(table(rec_df_10$hist_gra, rec_df_10$status_crr))
    # Set the multivariable columns and the tests
    dependent_os <- "Surv(rfs_upd, status_os)"
    dependent_dss <- "Surv(rfs_upd, status_dss)"
    dependent_crr <- "Surv(rfs_upd, status_crr)"
    
    or_cph <- rec_df_10 %>% 
        finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
        rename(!!or_label := label) %>% 
        rename(" " = levels) %>% 
        rename("  " = all)
    
    my_table <- rec_df_10 %>%
        # Summary table
        summary_factorlist(dependent_dss, explanatory, 
                           column = TRUE, fit_id = TRUE) %>%
        # CPH univariable
        ff_merge(
            explanatory %>%
                map(~ rec_df_10 %>%
                        coxph(as.formula(paste(dependent_dss, .x, sep = " ~ ")), 
                              data = .) %>%
                        fit2df(estimate_suffix = " CPH univariable")
                )%>%
                bind_rows()
        ) %>%
        # CPH multivariable
        ff_merge(
            rec_df_10 %>%
                coxphmulti(dependent_dss, explanatory) %>%
                fit2df(estimate_suffix = " CPH multivariable")
        ) %>%
        # Fine and Gray competing risks regression
        ff_merge(
            rec_df_10 %>%
                crrmulti(dependent_crr, explanatory) %>%
                fit2df(estimate_suffix = " competing risks multivariable")
        ) %>%
        select(-fit_id, -index) %>%
        dependent_label(rec_df_10, label)
    
    return(list(or_cph, my_table))
}

options(scipen = 0)
dir <- "~/PACC/TMA/Results/"
tables <- paste0(dir, "ManualScoreTables/")
plot_dir_surv <- paste0(dir, "SurvPlots/Survival/")
plot_dir_rec <- paste0(dir, "SurvPlots/Recurrence/")
patient_data <- read_tsv(paste0(dir,
                                "all_patient_unique_with_pacc_full_score_new.csv"),
                         na = "#NULL!") %>%
    mutate(PACC=if_else(PACC=="NA", NA, as.integer(PACC)),
           PACC_YN = if_else(PACC == 2, 1, 0),
           hist_gra=as.factor(hist_gra),
           recgrp=as.factor(recgrp),
           tumsize=as.numeric(tumsize),
           adjbeh = as.factor(adjbeh)
           
           ) %>% 
    filter(exkl == 0, tumor == 1) %>% 
    droplevels()

# Use 15 years censor for survival
# Year cut off to censor
censor_yr <- 15

## Survival ----
all_grades <- surv_hr(patient_df = patient_data,
                     label = "Survival")

write_tsv(all_grades[[1]], paste0(dir,"HazardRatio/os_with_grade.csv"))
write_tsv(all_grades[[2]], paste0(dir,"HazardRatio/surv_with_grade.csv"))

for (i in c(1,2,3)) {
    new_table <- surv_hr(patient_df = patient_data,
                         grade = i,
                         label = paste("Survival: Grade",i))
    
    write_tsv(new_table[[1]],
              paste0(dir,"HazardRatio/os_with_grade",i,".csv"))
    write_tsv(new_table[[2]],
              paste0(dir,"HazardRatio/surv_with_grade",i,".csv"))
}

## Overall and Local Recurrence ----
# Local recurrence
all_grades <- rec_hr(patient_data, grade = 0, 
                              label = "Local recurrence",
                              local = T)
gt(all_grades[[1]])

write_tsv(all_grades[[1]], paste0(dir,"HazardRatio/or_with_grade.csv"))
#write_tsv(all_grades[[2]], paste0(dir,"HazardRatio/lr_with_grade.csv"))
# Local recurrence for all grades separately
for (i in c(1,2,3)) {
    new_table <- rec_hr(patient_data, grade = i,
                              label = paste("Local recurrence: Grade", i),
                              local = TRUE)

    write_tsv(new_table[[1]], paste0(dir,"HazardRatio/or_with_grade",i,".csv"))
    #write_tsv(new_table[[2]], paste0(dir,"HazardRatio/lr_with_grade",i,".csv"))
}



## Distant recurrence table ----
# all_grades <- rec_hr(patient_data, grade = 0, 
#                               label = "Distant recurrence",
#                               local = F)
# write_tsv(all_grades[[2]], paste0(dir,"HazardRatio/dr_with_grade.csv"))
# for (i in c(1,2,3)) {
#     dist_table <- rec_hr(patient_data, grade = i, 
#                                  label = paste("Distant recurrence: Grade", i),
#                                  local = FALSE)
#     
#     write_tsv(dist_table[[2]], paste0(dir,"HazardRatio/dr_with_grade",i,".csv"))
# }
