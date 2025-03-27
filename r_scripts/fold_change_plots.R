library(tidyverse)
library(RColorBrewer)
library(ggbeeswarm)

setwd("/home/aboffelli/PACC/TMA/Results/")
##---- Human
human_score <- read_tsv("ManualScoreTables/ecc_fold_change.txt")
n_cells <- table(human_score$Pathologist)
gmm_table <- read_tsv("QuPathScoreTables/gmm_table.txt")
mean_cps <-gmm_table %>% 
    summarize(
        mean_cp1 = mean(mean_cp1),
        mean_cp2 = mean(mean_cp2)) %>% 
    mutate(foldchange= mean_cp1/mean_cp2)

all_cells <- read_tsv("QuPathScoreTables/all_cells_area.txt") %>% 
    mutate(Parent=str_remove(Parent,"_measurements.txt"),
           Image = str_remove_all(Image, "-|\\.tif")
           )



set.seed(1234)
all_cells <- all_cells %>% 
    slice_sample(n = 500) %>% 
    left_join(gmm_table, by = c("Parent", "Image")) %>% 
    mutate(FoldChange=Area/mean_cp2) 
    

fold_change_plot <- human_score %>%
    ggplot(aes(x=PaccArea, y=FoldChange, fill = Pathologist, 
               shape=Pathologist)) +
    geom_point(size=5) +
    geom_point(data=all_cells, aes(x=Area, y=FoldChange, 
                                      fill='#619CFF'),
               alpha=0.5, shape = 21, size=3) +
    geom_hline(yintercept = c(mean_cps$foldchange,1), linetype = c("dashed","dotdash"), 
               color = c("#A50026","#E08214"), size=2) +
    labs(x = "Nuclear Area", y = 'Area Fold Change', fill= 'Type', shape="Type", 
         title = "")+
    scale_shape_manual(values = c(21, 24),
                       labels = c(paste0("ECCs Pathologist 1 (n=", 
                                        n_cells["Pathologist 1"], ")"),
                                  paste0("ECCs Pathologist 2 (n=", 
                                               n_cells["Pathologist 2"], ")")
                                  )) +
    scale_fill_manual(values = c('#619CFF','#FFFF99', "#CAB2D6"))+
    scale_y_continuous(breaks = seq(0, 9, by = 1)) +
    guides(fill="none" , # Remove the fill legend
           shape = guide_legend(override.aes = list(fill = c('#FFFF99', "#CAB2D6")))) +  
    theme_classic(base_size = 30)+
    theme(legend.position = c(0.82,0.26),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)); test


ggsave("~/PACC/TMA/Results/pathologist_scatterplot.tif", fold_change_plot, width = 390, height=240, units = "mm")
