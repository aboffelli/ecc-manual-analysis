##---- Competing risk regression

# Suppress messages and warnings from the whole script to keep snakemake output clean
suppressMessages(suppressWarnings({
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

            # Remove variables with only one level, since not all grades have the same variables
            explanatory <- explanatory[
            vapply(explanatory, function(var) n_distinct(patient_df[[var]]) > 1, logical(1))
            ]                   
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
                
                recgrp = recode(as.character(recgrp),
                    "1" = "Luminal A",
                    "2" = "Luminal B",
                    "3" = "HER2",
                    "4" = "TN") %>%
                factor() %>%
                ff_label("Cancer subtype"),
                tumsize=ff_label(as.numeric(tumsize), "Tumor size"),
                PACC_YN=factor(PACC_YN) %>% 
                    fct_recode("-"= "0",
                            "+" = "1") %>% 
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

            # Remove variables with only one level, since not all grades have the same variables
            explanatory <- explanatory[
            vapply(explanatory, function(var) n_distinct(patient_df[[var]]) > 1, logical(1))
        ]
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


    # Access the snakemake variables
    input_file = snakemake@input[["patient_data"]]

    # Output files for survival
    out_surv = snakemake@output[["hr_surv"]]

    ## TODO: Fix the script to handle all grades separately
    out_surv_grade1 = snakemake@output[["hr_surv_grade_1"]]
    out_surv_grade2 = snakemake@output[["hr_surv_grade_2"]]
    out_surv_grade3 = snakemake@output[["hr_surv_grade_3"]]

    # Map the output files to a vector
    out_surv_grades <- c(out_surv_grade1, out_surv_grade2, out_surv_grade3)

    # Output files for recurrence
    out_rec = snakemake@output[["hr_rec"]]
    out_rec_grade1 = snakemake@output[["hr_rec_grade_1"]]
    out_rec_grade2 = snakemake@output[["hr_rec_grade_2"]]
    out_rec_grade3 = snakemake@output[["hr_rec_grade_3"]]
    # Map the output files to a vector
    out_rec_grades <- c(out_rec_grade1, out_rec_grade2, out_rec_grade3)

    options(scipen = 0)

    patient_data <- read_tsv(input_file,
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

    ## Survival ----
    # Survival with all grades
    all_grades <- surv_hr(patient_df = patient_data,
                        label = "Survival")

    write_tsv(all_grades[[2]], out_surv)

    # Survival for all grades separately
    for (i in c(1,2,3)) {
        new_table <- surv_hr(patient_df = patient_data,
                            grade = i,
                            label = paste("Survival: Grade",i))
        write_tsv(new_table[[2]], out_surv_grades[i])
    }

    ## Overall Recurrence ----
    # Recurrence with all grades
    all_grades <- rec_hr(patient_data, grade = 0, 
                                label = "Local recurrence",
                                local = T)

    write_tsv(all_grades[[1]], out_rec)

    # Recurrence for all grades separately
    for (i in c(1,2,3)) {
        new_table <- rec_hr(patient_data, grade = i,
                                label = paste("Local recurrence: Grade", i),
                                local = TRUE)

        write_tsv(new_table[[1]], out_rec_grades[i])
    }
}))