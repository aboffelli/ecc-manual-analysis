suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(ggbeeswarm)))


# Access the input files
ecc_fold_change <- snakemake@input[["ecc_fold_change"]]
gmm <- snakemake@input[["gmm_table"]]
all_cells_area <- snakemake@input[["all_cells_area"]]

# Access the output
output_file <- snakemake@output[[1]]

# Load the tables
human_score <- suppressMessages(suppressWarnings(
    read_tsv(ecc_fold_change)))
gmm_table <- suppressMessages(suppressWarnings(
    read_tsv(gmm)))

# Calculate the mean of each component for the reference lines
mean_cps <-gmm_table %>% 
    summarize(
        mean_cp1 = mean(mean_cp1),
        mean_cp2 = mean(mean_cp2)) %>% 
    mutate(foldchange= mean_cp1/mean_cp2)

# Load the tables containing all the detections from QuPath
all_cells <- suppressMessages(suppressWarnings(
    read_tsv(all_cells_area))) %>% 
    mutate(Parent=str_remove(Parent,"_measurements.txt"),
           Image = str_remove_all(Image, "-|\\.tif")
           )



set.seed(1234)  # Ensure reproducibility

# Since there are so many detections, we subset only 500 randomly from the data 
all_cells <- all_cells %>% 
    slice_sample(n = 500) %>% 
    left_join(gmm_table, by = c("Parent", "Image")) %>% # Join the gmm data
    mutate(FoldChange=Area/mean_cp2) # Calculate the fold change

## Plot ----

# Violin plot of the fold change by pathologist
fold_change_plot <- human_score %>% 
    ggplot(aes(x=Pathologist, y=FoldChange, fill=Pathologist))+
    geom_violin(trim=FALSE) +
    scale_fill_manual(values = c('#FFFF99', "#CAB2D6"))+  # Fill colors
    theme_classic()

# Scatter plot with fold change and area of the detections
# fold_change_plot <- human_score %>%
#     ggplot(aes(x=PaccArea, y=FoldChange, fill = Pathologist, 
#                shape=Pathologist)) +
#     geom_point(size=5) +    # Pathologist points
#     geom_point(data=all_cells, aes(x=Area, y=FoldChange, fill='#619CFF'),  
#                alpha=0.5, shape = 21, size=3) + # QuPath detection points smaller
#     geom_hline(yintercept = c(mean_cps$foldchange,1), linetype = c("dashed","dotdash"), 
#                color = c("#A50026","#E08214"), size=2) + # Reference lines
#     labs(x = "Nuclear Area", y = 'Area Fold Change', fill= 'Type', shape="Type", 
#          title = "")+   # Labels
#     scale_shape_manual(values = c(21, 24),  # Change shapes of points
#                        labels = c("ECCs Pathologist 1 (n=96)", 
#                                   "ECCs Pathologist 2 (n=134)")) +  # Legend labels
#     scale_fill_manual(values = c('#619CFF','#FFFF99', "#CAB2D6"))+  # Fill colors
#     scale_y_continuous(breaks = seq(0, 9, by = 1)) +    # Y axis breaks
#     guides(fill="none" , # Remove the fill legend (it breaks because of the all cells dots)
#            shape = guide_legend(override.aes = list(fill = c('#FFFF99', "#CAB2D6")))) +  # Add fill legend by shape
#     
#     # Increase the font and reposition the legend inside the plot
#     theme_classic(base_size = 30)+
#     theme(legend.position = c(0.82,0.26),
#           legend.title = element_text(size = 18),
#           legend.text = element_text(size = 18))

# Save the file
ggsave(output_file, fold_change_plot, width = 390, height=240, units = "mm")
