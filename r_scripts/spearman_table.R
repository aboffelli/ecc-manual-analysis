#---- Spearman correlation tests.
library(tidyverse)
library(ggpubr)

calculate_percent <- function(column) {
    # Function to calculate the percentage of patients by an specific category passed as an argument.
    
    group_summary <- patient_data %>%
        filter(!is.na({{ column }})) %>%  # remove NAs
        group_by({{ column }}, PACC) %>%  # Group by category and ECC score
        summarise(
            n = n(),  # Count of patients with this score in the group
            .groups = "drop") %>% 
        group_by({{ column }}) %>% # group again only by the category
        mutate(
            total_in_group = sum(n),  # Count the total amount of patients in that group.
            percentage = round((n / total_in_group) * 100, 0)  # calculate percentage
        ) %>%
        ungroup()
    
    return(group_summary)  
}

# Set up the input and output directories
dir <- "~/PACC/TMA/Results/"
tables <- paste0(dir, "ManualScoreTables/")
plot_dir_surv <- paste0(dir, "SurvPlots/OnlyGrade1and2/")
plot_dir_rec <- paste0(dir, "SurvPlots/OnlyGrade1and2/")

patient_data <- read_tsv(
    paste0(dir,
           "all_patients_with_TMA_PACC_HIF.csv")) %>%
    mutate(PACC=as.integer(PACC)) %>% 
    filter(PACC %in% c(0,1,2)) %>%  # Filter out ECC score NA
    droplevels() %>% 
    mutate(age_gap = if_else(age >= 50, ">=50", "<50"),  # Merge ages 
           tumor_size = if_else(tumsize <= 20, "<=20mm", ">20mm")  # Merge tumor sizes
           )  

## All ----

# Calculate the total amount of patients per ECC group.
stratified_summary <- patient_data %>%
    group_by(PACC) %>%
    summarise(
        n = n(),
        percentage = round(n() / nrow(patient_data) * 100,0)
    ) %>%
    ungroup()
print(stratified_summary)

## Age ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(age_gap)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$age, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Tumour size ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(tumor_size)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$tumsize, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Ki67 ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(ki67_pos)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$ki67_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=spearman_test$p.value)
print(rho)
print(p_value)

# Filter only NHG 2 patients for checking group correlation with Ki67 fraction
ki67_statistics <- patient_data %>% 
    filter(hist_gra == 2)

# fit the model to check normality
model <- lm(ki67_fra ~ PACC_YN, data = ki67_statistics)
plot(model)  # diagnostic plots
shapiro.test(residuals(model))  # Normality test

# The data is not normally distributed, so we use a non-parametric test.
# The test can be performed directly in ggplot using ggpubr.

# Plot of the data using violin plots
ki67_plot <- ki67_statistics %>% 
    ggplot(aes(x = as.factor(PACC_YN), y = ki67_fra)) +
    geom_violin(aes(fill=as.factor(PACC_YN))) + # Violin plot
    
    stat_compare_means(method = "wilcox.test", # Calculate the Wilcox Rank-Sum Test
                       label.x = 1.25,   # Adjust the position
                       size = 8 # Adjust font size
                       ) +
    
    scale_fill_manual(values = c("#B3CDE3", "#FBB4AE")) +   # Change the colors
    
    # Change labels
    scale_x_discrete(labels = c("ECC -", "ECC +")) +    
    labs(x = "ECC score", y = "Ki67 fraction",
         title = "Ki67 fraction by ECC group in NHG 2 patients") +
    
    theme_classic(base_size = 24) + # Put a clean theme and increase the font.
    
    theme(legend.position = "none") # since the y axis already has the labels we can get rid of the color legend.

# Save the tif
ggsave(paste0(dir, "ki67_grade2.tif"), ki67_plot, width = 10, height =10, units = "in")

## NHG ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(hist_gra)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$hist_gra, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## ER ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(er_pos)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$er_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## PR ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(pgr_pos)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$pgr_pos, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## HER2 ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(her2_com)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$her2_com, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## HIF ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(HIF1A_ACTIVATION)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$HIF1A_ACTIVATION, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=round(spearman_test$p.value, 3))
print(rho)
print(p_value)

## Grade markers - atypia ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(atypia)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$atypia, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=spearman_test$p.value)
print(rho)
print(p_value)

# Bar plot
patient_data %>% 
    ggplot(aes(x=PACC, fill=as.factor(atypia))) +
    geom_bar(position = "dodge") +
    labs(fill="Atypia", x = "ECC score") +
    theme_classic()

## Grade markers - mitosis ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(mitosis)
print(group_summary)

# Bar plot
patient_data %>% 
    ggplot(aes(x=PACC, fill=as.factor(mitosis))) +
    geom_bar(position = "dodge") +
    labs(fill="Mitosis", x = "ECC score") +
    theme_classic()

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$mitosis, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=spearman_test$p.value)
print(rho)
print(p_value)

## Grade markers - tubular formation ----

# Use the function to calculate the percentage
group_summary <- calculate_percent(tubuli)
print(group_summary)

# Perform the spearman correlation test
spearman_test <- cor.test(patient_data$tubuli, patient_data$PACC, method = "spearman")
rho <- round(spearman_test$estimate, 3)
p_value <- c(p.value=spearman_test$p.value)
print(rho)
print(p_value)

# Bar plot
patient_data %>% 
    ggplot(aes(x=PACC, fill=as.factor(tubuli))) +
    geom_bar(position = "dodge") +
    labs(fill="Tubular formation", x = "ECC score") +
    theme_classic()
