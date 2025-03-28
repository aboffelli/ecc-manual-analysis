## Gene expression data
library(tidyverse)
library(RColorBrewer)
library(asht)
options(scipen=0)


# Read the patient data, and set #NULL! to be recognized as NA.
data <- read_tsv("~/PACC/TMA/Results/all_patients_with_TMA_PACC_HIF_GEX.csv",
                 na = "#NULL!") %>% 
    #filter(hist_gra %in% c(1,2)) %>% 
    mutate(hist_gra=as.factor(hist_gra))    # Make a factor

# For the statistical test, we need the table in the long format.
long_exp <- data %>% 
    # Remove patients that do not have ECC score.
    filter(PACC_YN != "NA") %>%

    # Change the 0 and 1 to "-" and "+"
    mutate(PACC_YN = if_else(PACC_YN == 0, "-", "+")) %>%

    # Pivot the table ( only the genes columns )
    pivot_longer(cols = c("cdk1","epas1","hif1a","mki67"), names_to = "gene", 
                 values_to = "expression")


# Check for normality of the expression data
results_normality <- long_exp %>%
    group_by(PACC_YN, gene) %>%
    summarize(p_value = shapiro.test(expression)$p.value, .groups = 'drop')

# Since the data is skewed, we use mann_whitney.
results_mann_whitney <-  long_exp %>%
    group_by(gene) %>%

    # Get each value in one column (W, CI, and p-value)
    summarize(w_value= wmwTest(expression ~ PACC_YN)$statistic,
              CI=paste0(
                  round(wmwTest(expression ~ PACC_YN)$conf.int[[1]],2),
                  ", ",
                  round(wmwTest(expression ~ PACC_YN)$conf.int[[2]],2)),
              p_value=wmwTest(expression ~ PACC_YN)$p.value,
              .groups = 'drop')

## Plot ----
plots <- long_exp %>% 
    # Filter the gene you want (if you want a specific one)
    #filter(gene=="epas1") %>% 
    ggplot(aes(x=PACC_YN, y=expression)) +

    # Points for every patient.
    geom_jitter(aes(color = PACC_YN), size = 0.5, width = 0.35) +

    # Box for each ECC score
    geom_boxplot(aes(fill=PACC_YN), alpha=0.6) +

    # Change the fill colors (boxes)
    scale_fill_manual(values = c("-"='#9ECAE1',"+"='#F8766D'), 
                      name="", )+
    
    # Change the color colors (dots)
    scale_color_manual(values = c("-"='#9ECAE1',"+"='#F8766D'),
                       name="")+

    # Change theme and increase the text size
    theme_classic(base_size = 20)+

    # Remove legend, since the x-axis alread shows the group name.
    theme(legend.position = "none")+

    # Add the mann-whitney p-value on the top of the plot
    annotate("text", x=2, y= 1.5, label = paste0("p = ", round(results_mann_whitney$p_value, 3)), size=6) +

    # Change the labels (depends on the gene you filtered)
    labs(x="ECC Score", y="EPAS1 Expression")

    # If you want all the genes, uncomment the line below
    #facet_grid(~gene)

# Save a tif image.
ggsave("~/PACC/TMA/Manuscript/gex_boxplot.tif", plots, width=110, height=110, units = 'mm')
