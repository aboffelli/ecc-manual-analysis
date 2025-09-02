##---- Kaplan meier curves PACC / no-PACC

suppressMessages(suppressWarnings({
    library(tidyverse)
    library(survival)
    library(survminer)
    library(ggsurvfit)
    library(gridExtra)
    library(tidycmprsk)
    library(gtsummary)
    library(RColorBrewer)

    # Suppress empty Rplots.pdf
    if (length(dev.list()) == 0) {
    pdf(NULL)
    }

    # options(scipen = 0) # Avoid scientific notation for decimals.

    # Fuction to generate the survival plot.
    survival_km <- function(patient_df, grade=0, label) {
        
        if (!grade %in% c(1,2,3)) {
            # If grade is not specified, just set a vector with all grades
            grade <- c(1,2,3, NA)
            label <- paste0(label,": All NHG")
        } else {
            # If the grade is specified, add with grade is being used in the label.
            label <- paste0(label,": NHG ", grade)   
        }
        
        censor_yr <- 15 # 15 years time limit for survival.
        
        surv_df <- patient_df %>% 
            mutate(hist_gra=as.integer(hist_gra),
                PACC_YN = if_else(PACC != 2, 0, 1),  # Make the ECC score binary
                status_os = if_else(bcd_cmp_ == 1, 1, 0)  # Only death by breast cancer
                ) %>% 

            filter(tumor == 1, hist_gra %in% grade) %>% # Filter only primary tumor and only the grades selected.
            droplevels()
        
        # Censoring the data in 15 years
        surv_df_15 <- surv_df %>%
            mutate(status_os = if_else(surv_dec > censor_yr, 0, status_os), # patients that died after 15 years should be marked as alive (0)
                surv_dec = if_else(surv_dec > censor_yr, censor_yr, surv_dec))   # if the follow up time is higher than 15 years, change to 15
        
        # Fit model
        fit <- survfit2(
            Surv(surv_dec, status_os) ~ PACC_YN,  # Fit the model using the ECC score
            data = surv_df_15)
        #summary(fit)
        
        # Statistics Rho=0 for log-rank, Rho=1 for Peto&Peto mod.
        logrank_test <- survdiff(Surv(surv_dec, status_os) ~ PACC_YN, 
                                data = surv_df_15, rho=0)

        # Create the plot
        km_plot <- fit %>%
            ggsurvfit() +

            # Change the labels
            labs(title = label,
                x = 'Time (years)') +

            # Add censor marks in the line
            add_censor_mark() +

            # Add the table under the plot
            add_risktable() +

            # Change the color and label for the lines
            scale_color_manual(values = c('#619CFF','#F8766D'),
                            labels = c('ECC-', 'ECC+')) +
            scale_fill_manual(values = c('#619CFF','#F8766D'),
                            labels = c('ECC-', 'ECC+'))+
            
            # Keep the y axis from 0 to 1 to be consistent between plots
            ylim(c(0,1)) +

            # Add the log rank value in the plot, under the lines to the left
            annotate("text", x = 2.5, y = 0.65,
                    label = paste0("Log-rank p = ", round(logrank_test$pvalue, 4))) +
            
            # Add a clean theme and increase the text font
            theme_classic() +
            theme(text = element_text(size = 16))

        return(km_plot) # return only the plot
    }

    # Fuction to generate the recurrence free survival plot.
    recurrence_km <- function(patient_df, grade=0, 
                            label) {

        if (!grade %in% c(1,2,3)) {
                grade <- c(1,2,3, NA)
                label <- paste0(label,": All NHG")
        } else {
                label <- paste0(label,": NHG ", grade)   
            }

        censor_yr <- 10     # Censor in 10 years for RFS

        # Prepare the data
        rec_df <- patient_df %>%

            # Remove the type 5 which are unknown recurrence
            filter(tumor == 1, ev1typ20 != 5) %>%

            mutate(hist_gra=as.integer(hist_gra),
                or_rec = if_else(ev1typ20 != 0, 1, 0), # Make the recurrence binary, any type of recurrence is 1
                PACC_YN = if_else(PACC != 2, 0, 2),  # Make the ECC score binary
            ) %>%
            
            # Filter the grade 
            filter(hist_gra %in% grade) %>% 
            droplevels()

        # Censor all the data in 10 years instead of full time.
        rec_df_10 <- rec_df %>%
            mutate(or_rec = if_else(rfs_ny > censor_yr, 0, or_rec), # If RSF follow up is higher than 10 years, change to no recurrence
                rfs_upd = if_else(rfs_ny > censor_yr, censor_yr, rfs_ny))    # If RSF follow up is higher than 10 years, change to 10 years

        # Fit the survival model
        fit_rec <- survfit2(
            Surv(rfs_upd, or_rec) ~ PACC_YN, # Fit the model using the ECC score
            data = rec_df_10)

        # Calculate the log-rank value
        logrank_rec <- survdiff(Surv(rfs_upd, or_rec) ~ PACC_YN,
                                data = rec_df_10, rho=0)

        # Create the plot
        km_plot <- fit_rec %>%
            ggsurvfit() +
            
            # Change the labels
            labs(title = label,
                x = 'Time (years)',
                y = 'Recurrence-free Survival Probability') +

            # Add censor marks in the line
            add_censor_mark() +

            # Add the table under the plot
            # TODO: Make the risktable work with the multiple plots pdf
            add_risktable() +

            # Change the color and label for the lines
            scale_color_manual(values = c('#619CFF','#F8766D'),
                            labels = c('ECC-', 'ECC+')) +
            scale_fill_manual(values = c('#619CFF','#F8766D'),
                            labels = c('ECC-', 'ECC+'))+

            # Keep the y axis from 0 to 1 to be consistent between plots
            ylim(c(0,1)) +

            # Add the log rank value in the plot, under the lines to the left
            annotate("text", x = 2.5, y = 0.65,
                    label = paste0("Log-rank p = ", round(logrank_rec$pvalue, 4))) +

            # Add a clean theme and increase the text font      
            theme_classic() +
            theme(text = element_text(size = 16))

        return(km_plot) # return just the plot
    }

    # Access the input files
    patient_data <- snakemake@input[["patient_data"]]

    # Access output files
    surv_dist <- snakemake@output[["surv_dist"]]
    out_surv <- snakemake@output[["surv_curves"]]
    out_rec <- snakemake@output[["rec_curves"]]

    # Load the patient data
    patient_data <- read_tsv(patient_data, na="#NULL!") %>%
        mutate(PACC=as.integer(PACC)) %>% 
        
        filter(exkl == 0,  # remove patients that were excluded from the study
            PACC %in% c(0,1,2),  # remove ECC score NA
            tumor==1) %>%    # Keep only primary tumor data

        droplevels()

    ## Data exploration plots ----

    # Survival distribution
    plot1 <- patient_data %>% 
        ggplot(aes(x=surv_dec)) +   # Follow up time
        geom_histogram(aes(fill=as.factor(bcd_cmp_)),   # Death variable
                    col='black', alpha=0.5) +
        
        # Change the labels, and remove the legend title since it's not necessary
        labs(title="Distribution of survival time", x="Survival Time (years)", 
            y='Count', fill="") +

        # Change the colors and labels for the death types
        scale_fill_manual(
            values = c("#1B9E77", "#D95F02", "#7570B3"),
            labels = c("Alive",     # Value 0 - Alive
                    "Death by breast cancer",    # Value 1 - Breast cancer death
                    'Death by other means')  # Value 2 - Other deaths
            ) +
        
        # Add a clean theme and increase the font size
        theme_classic(base_size = 20)

    # Save a tif file
    ggsave(surv_dist, plot1, width = 9, 
        height = 9, 
        units = "in")
        

    ## Survival plots ----

    # Create the loop to generate one plot for each NHG. 0 represents all grades.
    plots <- list()
    for (i in c(0,1,2,3)) {
        surv_plot <- survival_km(patient_df = patient_data, 
                                label = "Kaplan-Meier Survival Curve", grade = i)
        plots[[length(plots) + 1]] <- surv_plot
    }
    ggsave(out_surv, marrangeGrob(plots, nrow = 1, ncol = 1))

    ## Recurrence plots ----

    # Create the loop to generate one plot for each NHG. 0 represents all grades. 
    plots <- list() # List to store the plots
    for (i in c(0,1,2,3)) {
        # Use the function to generate the plots
        rec_plot <- recurrence_km(patient_df = patient_data, 
                                    label="Kaplan-Meier RFS Curve", grade=i)
        plots[[length(plots) + 1]] <- rec_plot
    }
    ggsave(out_rec, marrangeGrob(plots, nrow = 1, ncol = 1))

})) # End of suppress messages and warnings