library(tidyverse)
library(gridExtra)
library(grid)
library(mclust)
library(gtools)

setwd(dir = "~/PACC/TMA/Results/QuPathScoreTables/")

plot_gmm_histogram <- function(data, image_name, filename) {
    filename = str_remove(filename, "_measurements.txt")
    print(paste(filename, image_name))
    gmm <- Mclust(data$`Area µm^2`, G = 2)
    

    means <- gmm$parameters$mean
    variances <- gmm$parameters$variance$sigmasq
    proportions <- gmm$parameters$pro
    
    # Create a sequence of points across the data range to plot the GMM line
    x_seq <- seq(min(data$`Area µm^2`), max(data$`Area µm^2`), 
                 length.out = 1000)
    
    # Calculate the GMM density for each component
    gmm_density <- rep(0, length(x_seq))

    # Loop through each component of the GMM
    for (k in 1:gmm$G) {
        gmm_density <- gmm_density + proportions[k] * 
            dnorm(x_seq, mean = means[k], 
                  sd = sqrt(variances[k]))
    }
    
    # Prepare the data frame for plotting
    density_data <- data.frame(x = x_seq, y = gmm_density)
    

    p <- data %>% 
        ggplot(aes(`Area µm^2`)) +
        geom_histogram(aes(y = after_stat(density)), bins = 300, 
                       fill = "#619CFF", color = "#08519C", alpha=0.5) +
        geom_line(data = density_data, aes(x = x, y = y), color = "black",
                  linewidth = 1.2) +
        geom_vline(xintercept = means, linetype = c("dashed", "dotdash"), ,
                   color = c("#A50026", "#E08214"), size = 2) +
        ggtitle(image_name) +
        labs(x="Nuclear Area", y="Density") +
        theme_classic(base_size = 20)+
        xlim(0,1000)
    
    print(p)

    # Return summary statistics as a tibble
    tibble(
        Parent = filename,
        Image = image_name,
        mean_cp1 = means[1],
        mean_cp2 = means[2],
        variation_cp1 = variances[1],
        variation_cp2 = variances[2],
        proportion_cp1 = proportions[1],
        proportion_cp2 = proportions[2]
    )

}

filenames <- c("1A_EPCAM_measurements.txt",
               "2A_EPCAM_measurements.txt",
               "15A_EPCAM_measurements.txt",
               "20A_EPCAM_measurements.txt")

results <- filenames %>%
    map_dfr(~ {
        data <- read_tsv(.x) %>% 
            mutate(Image=as.character(str_replace_all(Image, "-|\\.tif$", "")))
        
        # Create a PDF for the current file
        pdf_filename <- str_replace(.x, "_measurements.txt", "_histograms.pdf")
        pdf(pdf_filename,  width = 390 / 25.4, height = 240 / 25.4)
        
        # Process each Image within the file and plot
        image_results <- data %>%
            group_by(Image) %>%
            do(plot_gmm_histogram(., unique(.$Image), .x))
        dev.off()
        image_results
    })

results %>%
    write_tsv("gmm_table.txt")
