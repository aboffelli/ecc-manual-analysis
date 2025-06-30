source("renv/activate.R")
## Gene expression data
library(tidyverse)
library(RColorBrewer)
library(asht)
options(scipen=0)


dir <- "~/PACC/TMA/Results/"
data <- read_tsv(paste0(dir,
                        "all_patients_with_TMA_PACC_HIF_GEX.csv"),
                 na = "#NULL!") %>% 
    #filter(hist_gra %in% c(1,2)) %>% 
    mutate(hist_gra=as.factor(hist_gra))
long_exp <- data %>% 
    filter(PACC_YN!="NA") %>% 
    pivot_longer(cols = c("cdk1","epas1","hif1a","mki67"), names_to = "gene", values_to = "expression")


# Check for normality of the expression data
results_normality <- long_exp %>%
    group_by(PACC_YN, gene) %>%
    summarize(p_value = shapiro.test(expression)$p.value, .groups = 'drop')


results_mann_whitney <-  long_exp %>%
    group_by(gene) %>%
    summarize(w_value= wmwTest(expression ~ PACC_YN)$statistic,
              CI=paste0(
                  round(wmwTest(expression ~ PACC_YN)$conf.int[[1]],2),
                  ", ",
                  round(wmwTest(expression ~ PACC_YN)$conf.int[[2]],2)),
              p_value=wmwTest(expression ~ PACC_YN)$p.value,
              .groups = 'drop') %>% 
    filter(gene=="hif1a")

plots <- long_exp %>% 
    filter(gene=="hif1a") %>% 
    ggplot(aes(x=PACC_YN, y=expression)) +
    geom_jitter(aes(color = PACC_YN), size = 0.5, width = 0.35) +
    geom_boxplot(aes(fill=PACC_YN), alpha=0.6) +
    scale_fill_manual(values = c("0"='#9ECAE1',"2"='#F8766D'), 
                      name="", )+
    scale_color_manual(values = c("0"='#9ECAE1',"2"='#F8766D'),
                       name="")+
    theme_classic(base_size = 20)+
    theme(legend.position = "none")+
    labs(x="ECC Score", y="HIF1A Expression"); plots
#facet_grid(~gene); plots

ggsave("~/PACC/TMA/Manuscript/gex_boxplot.tif", plots, width=110, height=110, units = 'mm')
