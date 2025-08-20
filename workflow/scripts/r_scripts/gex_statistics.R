## Gene expression data

# Suppress messages and warnings from the whole script to keep snakemake output clean
suppressMessages(suppressWarnings({
    library(tidyverse)
    library(RColorBrewer)
    if (!requireNamespace("asht", quietly = TRUE)) {
        install.packages("asht", repos = "https://cloud.r-project.org")
        }
    library(asht)
    options(scipen=0)


    # Access the input files
    input_file <- snakemake@input[[1]]

    # Access output
    output_file <- snakemake@output[[1]]


    # Read the patient data, and set #NULL! to be recognized as NA.
    data <- suppressMessages(suppressWarnings(read_tsv(input_file,
                    na = "#NULL!"))) %>% 
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
    gene_to_plot <- snakemake@params[["gene"]]
    y_label <- snakemake@params[["y_label"]]

    plots <- long_exp %>% 
        # Filter the gene you want (if you want a specific one)
        filter(gene==gene_to_plot) %>% 
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
        annotate("text", x=1.5, y= 2, label = paste0("p = ", round(results_mann_whitney$p_value[results_mann_whitney$gene==gene_to_plot], 3)), size=6) +

        # Change the labels (depends on the gene you filtered)
        labs(x = "ECC Score", y = y_label)

        # If you want all the genes, uncomment the line below
        #facet_grid(~gene)

    # Save a tif image.
    ggsave(output_file, plots, width=110, height=110, units = 'mm')

}))
# End of suppressMessages