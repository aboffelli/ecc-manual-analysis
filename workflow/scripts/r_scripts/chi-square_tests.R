source("renv/activate.R")
## Chi-square tests
library(tidyverse)
library(car)    # For ANOVA
library(emmeans)
library(stats)
library(chisq.posthoc.test)
library(knitr)
library(gt)
library(RColorBrewer)
library(dunn.test)
library(mice)
library(VIM)

options(scipen = 0)
dir <- "~/TMA/Results/"
tables <- paste0(dir, "ManualScoreTables/")
plot_dir_surv <- paste0(dir, "SurvPlots/Survival/")
plot_dir_rec <- paste0(dir, "SurvPlots/Recurrence/")
data <- read_tsv(paste0(dir,
                                "all_patient_unique_with_pacc_full_score.csv"),
                 na = "#NULL!") %>%
    mutate(PACC= if_else(PACC == "NA", NA, as.integer(PACC)),
           PACC=as.factor(PACC),
           hist_gra=as.factor(hist_gra),
           recgrp_h=as.factor(recgrp_h),
           tumsize=as.numeric(tumsize),
           adjbeh = as.factor(adjbeh),
           PACC_YN = if_else(PACC != 0, 1, 0)) %>% 
    filter(exkl == 0) %>% 
    droplevels()


# Step 2: Create a contingency table (countas.integer()# Step 2: Create a contingency table (count matrix)
contingency_table <- table(data$PACC_YN, data$hist_gra)
print("Contingency Table:")
print(contingency_table)

# Step 3: Calculate the expected matrix
expected_matrix <- chisq.test(contingency_table)$expected
print("Expected Matrix:")
print(expected_matrix)

# Step 4: Perform chi-square test
result <- chisq.test(contingency_table)
print("Chi-square Test Result:")
print(result)

posthoc_result <- chisq.posthoc.test(contingency_table, method = "bonferroni",
                                     round = 4) %>% 
    rename(`1`="1", `2`="2", `3`="3")

print(posthoc_result)

posthoc_result %>% 
    filter(Value=="Residuals") %>% 
    pivot_longer(cols = c("1","2","3"), 
                 names_to = "Hist_gra", values_to = "Residuals") %>% 
    ggplot(aes(x=Dimension, y=Hist_gra, fill=Residuals)) +
    scale_fill_distiller(palette = "RdBu") +
    geom_tile() +
    labs(x="PACC score", y="Histology grade") +
    theme_minimal()

# Convert the table to a graphical object
# my_table <- posthoc_result %>%
#     rename("PACC score"=Dimension) %>% 
#     write.table(file = "~/Desktop/table_data.csv", row.names = FALSE, quote = F,
#                 sep = '\t')

# Calculate the time to death after local recurrence
time_to_death <- data %>% 
    filter(ev1typ20 == 1, bcd_cmp_ == 1, tumor == 1) %>% 
    select(ev1typ20, bcd_cmp_, PACC, hist_gra, rfs_upd, surv_dec) %>% 
    mutate(time_to_death=surv_dec-rfs_upd) %>% 
    drop_na()

grade_plot <- time_to_death %>% 
    ggplot(aes(y=time_to_death, x=hist_gra, fill=hist_gra)) +
    geom_boxplot() + 
    geom_jitter(width=0.2) +
    labs(x='Grade', y="Time",
         title = "Grade vs Time to death after local recurrence") +
    guides(fill = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic(); grade_plot
pacc_plot <- time_to_death %>% 
    ggplot(aes(y=time_to_death, x=PACC, fill=PACC)) +
    geom_boxplot()+
    geom_jitter(width=0.2) +
    labs(x='PACC score', y="Time",
         title = "PACC vs Time to death after local recurrence") +
    guides(fill = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic(); pacc_plot

# Calculate the time to death after distant recurrence
time_to_death_dist <- data %>% 
    filter(ev1typ20 == 4, bcd_cmp_ == 1, tumor == 1) %>% 
    select(ev1typ20, bcd_cmp_, PACC, hist_gra, rfs_upd, surv_dec) %>% 
    mutate(time_to_death=surv_dec-rfs_upd) %>% 
    drop_na()

grade_plot_dist <- time_to_death_dist %>% 
    ggplot(aes(y=time_to_death, x=hist_gra, fill=hist_gra)) +
    geom_boxplot() + 
    geom_jitter(width=0.2) +
    labs(x='Grade', y="Time",
         title = "Grade vs Time to death after distant recurrence") +
    guides(fill = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic(); grade_plot_dist
pacc_plot_dist <- time_to_death_dist %>% 
    ggplot(aes(y=time_to_death, x=PACC, fill=PACC)) +
    geom_boxplot()+
    geom_jitter(width=0.2) +
    labs(x='PACC score', y="Time",
         title = "PACC vs Time to death after distant recurrence") +
    guides(fill = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic(); pacc_plot_dist

# Save the plots
# ggsave(filename = "~/Desktop/grade_time_local.png", plot = grade_plot,
#        width=13, units = 'cm')
# ggsave(filename = "~/Desktop/pacc_time_local.png", plot = pacc_plot,
#        width=13, units = 'cm')
# ggsave(filename = "~/Desktop/grade_time_dist.png", plot = grade_plot_dist,
#        width=13, units = 'cm')
# ggsave(filename = "~/Desktop/pacc_time_dist.png", plot = pacc_plot_dist,
#        width=13, units = 'cm')

anova_result <- aov(time_to_death ~ PACC, data = time_to_death_dist)

# Display ANOVA table
summary(anova_result)

# Perform post-hoc tests (e.g., Tukey's HSD)
posthoc_results <- emmeans(anova_result, specs = pairwise ~ PACC, 
                           adjust = "tukey")

# Display post-hoc test results
pairs(posthoc_results)

# 1. Normality of residuals
residuals <- residuals(anova_result)
qqnorm(residuals)
qqline(residuals)
# You can also use Shapiro-Wilk test
shapiro.test(residuals)

# 2. Homogeneity of variances
# Plot residuals against predicted values or grouping variable
plot(anova_result, which = 1)  # Residuals vs. Fitted
# Use Levene's test
leveneTest(anova_result)

# 3. Independence of Errors
# Plot residuals against other variables (if applicable)
plot(anova_result, which = 2)  # Scale-Location plot
# Check for autocorrelation
acf(residuals)

# Since the residuals are not normally distributed, we use Kruskal-Wallis test.
kruskal_result_gra <- kruskal.test(time_to_death ~ hist_gra, 
                               data = time_to_death)
print(kruskal_result_gra)
kruskal_result_pacc <- kruskal.test(time_to_death ~ PACC, 
                                    data = time_to_death)
print(kruskal_result_pacc)

# Kruscal-wallis distant recurrence
kruskal_result_gra_dist <- kruskal.test(time_to_death ~ hist_gra, 
                                   data = time_to_death_dist)
print(kruskal_result_gra_dist)
kruskal_result_pacc_dist <- kruskal.test(time_to_death ~ PACC, 
                                         data = time_to_death_dist)
print(kruskal_result_pacc_dist)

# PostHoc Dunn's test for distant rec PACC (the only significant)
dunn.test(time_to_death_dist$time_to_death, time_to_death_dist$PACC, 
          method = "bonferroni")

# Missing data
time_to_death <- data %>% 
    filter(ev1typ20 == 1, bcd_cmp_ == 1, tumor==1) %>% 
    select(ev1typ20, bcd_cmp_, PACC, hist_gra, rfs_upd, surv_dec) %>% 
    mutate(time_to_death=surv_dec-rfs_upd)

aggr(time_to_death, prop=T,cex.axis=0.6)

mi_table <- time_to_death %>% 
    select(PACC, hist_gra, time_to_death)

dataimputed <- mice(mi_table, m=20, maxit=100, seed=1142, print=FALSE)
dataimputed$imp$PACC
plot(dataimputed)

kruskal_results <- with(dataimputed, kruskal.test(time_to_death ~ PACC))
