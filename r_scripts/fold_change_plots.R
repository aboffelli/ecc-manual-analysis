library(tidyverse)
library(RColorBrewer)
library(ggbeeswarm)

##---- Human
human_score <- read_tsv("/home/aboffelli/PACC/TMA/Results/ManualScoreTables/new_fixed_manual_scoring_sizes.csv") %>% 
    mutate(FoldChange= (PaccArea/CoreMean))
all_cells <- read_tsv("/media/aboffelli/ArthurExternalHD/BackUp/TMA/Data/QuPathImages/EPCAM/merged_measurements.txt")

#human_score$TMA <- sub("2", "15", human_score$TMA)

mean <- mean(human_score$FoldChange)
median <- median(human_score$FoldChange)
std_dev <- sd(human_score$FoldChange)

core_means <- human_score %>% 
    select(TMA,Core,CoreMean,BottomCI,UpperCI)

set.seed(1234)
merged_table <- merge(all_cells, core_means, by = c("TMA", "Core"))
merged_table <- merged_table %>% 
    mutate(InCI=ifelse(merged_table$`Area µm^2` >= merged_table$BottomCI &
                           merged_table$`Area µm^2` <= merged_table$UpperCI, 
                       1,0),
           FoldChange = (`Area µm^2`/CoreMean)) %>% 
    #filter(TMA=="2A") %>% 
    slice_sample(n = 500)
unique_foldchange <- as.data.frame(unique(round(merged_table$FoldChange, 2)))
colnames(unique_foldchange) <- "FoldChange" 

test <- human_score %>%
    ggplot(aes(x=PaccArea, y=FoldChange, fill = Pathologist, 
               shape=Pathologist)) +
    geom_point(size=5) +
    geom_point(data=merged_table, aes(x=`Area µm^2`, y=FoldChange, 
                                      fill='#619CFF'),
               alpha=0.5, shape = 21, size=3) +
    geom_hline(yintercept = c(0.5,1), linetype = c("dashed","dotdash"), 
               color = c("#A50026","#E08214"), size=2) +
    labs(x = "Nuclear Area", y = 'Area Fold Change', fill= 'Type', shape="Type", 
         title = "")+
    scale_shape_manual(values = c(24, 21),
                       labels = c("ECCs Pathologist 1 (n=96)", 
                                  "ECCs Pathologist 2 (n=134)")) +
    scale_fill_manual(values = c('#619CFF','#FFFF99', "#CAB2D6"))+
    guides(fill="none" , # Fix legend appearance
           shape = guide_legend(override.aes = list(fill = c('#FFFF99', "#CAB2D6")))) +  # Consistent shape and fill
    theme_classic(base_size = 30)+
    theme(legend.position = c(0.82,0.2),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)); test
    #facet_wrap(~ Core)#facet_wrap(~ Core)#facet_wrap(~ Core)

ggsave("~/PACC/TMA/Manuscript/pathologist_scatterplot.tif", test, width = 390, height=240, units = "mm")
# count_groups <- human_score %>% 
#     group_by(scale) %>% 
#     mutate(count = n())

# human_score %>% 
#     ggplot(aes(x=TMA, y=Means, col=TMA)) +
#     geom_boxplot() +
#     geom_point(position="jitter") +
#     theme_classic()

count_groups %>% 
    ggplot(aes(x=PACC_fc, y=count)) +
    geom_point(aes(col=TMA),position = "jitter") +
    #geom_density(aes(y=..count..), size = 1) +
    theme_classic() +
    scale_color_brewer(palette = "Paired") +
    geom_vline(xintercept = mean, linetype = "dashed", 
               color = "brown", size = 1, alpha=0.7) +
    geom_vline(xintercept = mean + std_dev, linetype = "dashed", 
               color = "sandybrown", size = 1, alpha=0.7) +
    geom_vline(xintercept = mean - std_dev, linetype = "dashed", 
               color = "sandybrown", size = 1, alpha=0.7) +
    #geom_vline(xintercept = median, linetype = "dotdash", 
     #          color = "brown", size = 1, alpha=0.7) +
    labs(x="Area fold change", y="Count", 
         title = "PACCs area fold change - Pathologist 2") +
    geom_text(aes(label = sprintf("Mean: %.2f\n N: %d", round(mean,2), 
                                  nrow(count_groups))),
              x = mean, y = 0, vjust = 0, hjust = -0.2, color = "brown", size = 4) +
    geom_text(aes(label=sprintf("SD: %.2f", mean + std_dev)), x = mean + std_dev,
              y = 0, vjust=0, hjust=-0.2, color="sandybrown", size=3) +
    geom_text(aes(label=sprintf("SD: %.2f", mean - std_dev)), x = mean - std_dev,
              y = 0, vjust=0, hjust=-0.2, color="sandybrown", size=3) +
    scale_x_continuous(limits = c(2, max(as.integer(count_groups$PACC_fc)) + 1)) +
    theme(text = element_text(size = 24))

##---- QuPath
ai_score <- read_tsv("~/qp_score") %>% 
    mutate(TMA=as.factor(TMA)) %>% 
    mutate(scale=as.factor(round(PACCS_fc)))

mean <- mean(ai_score$PACCS_fc)
median <- median(ai_score$PACCS_fc)
std_dev <- sd(ai_score$PACCS_fc)

count_groups <- ai_score %>% 
    group_by(scale) %>% 
    mutate(count = n())

# ai_score %>% 
#     ggplot(aes(x=TMA, y=Means, col=TMA)) +
#     geom_boxplot() +
#     geom_point(position="jitter") +
#     theme_classic()

count_groups %>% 
    ggplot(aes(x=PACCS_fc, y=count)) +
    geom_point(aes(col=TMA),position = "jitter") +
    #geom_density(aes(y=..count..), size = 1) +
    theme_classic() +
    scale_color_brewer(palette = "Paired") +
    geom_vline(xintercept = mean, linetype = "dashed", 
               color = "brown", size = 1, alpha=0.7) +
    geom_vline(xintercept = mean + std_dev, linetype = "dashed", 
               color = "sandybrown", size = 1, alpha=0.7) +
    geom_vline(xintercept = mean - std_dev, linetype = "dashed", 
               color = "sandybrown", size = 1, alpha=0.7) +
    # geom_vline(xintercept = median, linetype = "dotdash", 
    #            color = "brown", size = 1, alpha=0.7) +
    labs(x="Mean area fold change", y="Count", 
         title = "PACCs mean area fold change - QuPath") +
    geom_text(aes(label = sprintf("Mean: %.2f\n N: %d", round(mean,2), nrow(count_groups))),x = mean, y = 0, vjust = 0, hjust = -0.2, 
              color = "brown", size = 4) +
    geom_text(aes(label=sprintf("SD: %.2f", mean + std_dev)), x = mean + std_dev,
              y = 0, vjust=0, hjust=-0.2, color="sandybrown", size=3) +
    geom_text(aes(label=sprintf("SD: %.2f", mean - std_dev)), x = mean - std_dev,
              y = 0, vjust=0, hjust=-0.2, color="sandybrown", size=3) +
    scale_x_continuous(limits = c(2, 10)) +
    theme(text = element_text(size = 24))
    
## Log2 histogram

human_score_area <- human_score %>% 
    mutate(log2scale=log2(PaccArea),
           Type="PACCs") %>% 
    select(PaccArea, log2scale, Type) %>% 
    rename(Area=PaccArea)

set.seed(1689)
all_cells_area <- all_cells %>% 
    mutate(log2scale=log2(`Area µm^2`),
           Type="2n cells") %>% 
    select(`Area µm^2`, log2scale, Type) %>% 
    rename(Area=`Area µm^2`) %>% 
    slice_sample(n = 2000)

joined <- bind_rows(human_score_area, all_cells_area)
plot<- joined %>% 
    ggplot(aes(x=Area, fill = Type)) +
    geom_histogram(binwidth = 5, color = "black") +
    labs(x = "Area", title = "Histogram of Area with Log2 Transformed Counts") +
    theme_classic()

# Extract the counts from the histogram
hist_counts <- ggplot_build(plot)$data[[1]]

# Add a small constant to avoid log(0) issues and apply log2 transformation
hist_counts <- hist_counts %>%
    mutate(log2_count = log2(count + 1))

# Print the transformed data
print(hist_counts)

ggplot(hist_counts, aes(x = x, y = log2_count, fill=fill)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#D6604D","#4393C3"),
                      name = "",
                      labels=c("PACCs","2n cells"))+
    labs(x = "Area", y = "Log2(Count + 1)", 
         title = "Histogram of Nuclear Area with Log2 Transformed Counts") +
    theme_minimal()
