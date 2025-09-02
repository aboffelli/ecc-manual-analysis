library(tidyverse)
library(gridExtra)
library(grid)
library(mclust)
library(gtools)

setwd(dir = "~/PACC/TMA/Results/QuPathScoreTables/")

plot_gmm_histogram <- function(data, image_name, filename) {
    filename = str_remove(filename, "_measurements.txt") # get only the TMA number.
    print(paste(filename, image_name))  # Print just as a progress indicator
    
    gmm <- Mclust(data$`Area µm^2`, G = 2)  # calculate the Gaussian mixture
    
    # Calculate the percentiles for the Area.
    quartiles <- quantile(data$`Area µm^2`, probs = c(0.25, 0.50, 0.75, 0.99))
    q1 <- quartiles[1]  # 25th percentile
    median <- quartiles[2]  # 50th percentile (median)
    q3 <- quartiles[3]  # 75th percentile
    cutoff_99 <- quartiles[4]  # 99th percentile
    
    # Extract the means, variance, and proportion from the gm model.
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
        
        # Histogram for the area
        geom_histogram(aes(y = after_stat(density)), bins = 300, 
                       fill = "#619CFF", color = "#08519C", alpha=0.5) +
        
        # Line for the GMM fit
        geom_line(data = density_data, aes(x = x, y = y), color = "black",
                  linewidth = 1.2) +
        
        ## Add the quartiles lines (maybe not necessary)
        # geom_vline(xintercept = c(q1, median, q3), linetype = "dotted", 
        #            color = "darkgray", linewidth = 1) +  # Quartiles
        
        # Add the 99th percentile cutoff
        geom_vline(xintercept = cutoff_99, linetype = "dashed", 
                   color = "blue", linewidth = 1) +
        
        # Add the GMM mean lines for 1st and 2nd components.
        geom_vline(xintercept = means, linetype = c("dashed", "dotdash"), ,
                   color = c("#A50026", "#E08214"), size = 1) +
        
        ## Annotate the quartile lines.
        # annotate("text", x = c(q1, median, q3), y = 0.025, 
        #          label = c("Q1", "Median", "Q3"), angle = 90, 
        #          vjust = -0.5, size = 4, color = "darkgray") +
        
        # Annotate the 99th percentile line
        annotate("text", x = cutoff_99, y = 0.028, label = "99%", angle = 90, 
                 vjust = -0.5, size = 4, color = "blue") +
        
        # Annotate the GMM component lines
        annotate("text", x = means, y = 0.026, 
                 label = c("1st comp mean", "2nd comp mean"), angle = 90, 
                 vjust = -0.5, size = 4, color = c("#A50026", "#E08214")) +
        
        # Edit the title and labels
        labs(title = image_name,
             x = "Nuclear Area (µm²)",
             y = "Density") +
        
        # clean theme
        theme_classic(base_size = 20)+
        
        # Increase the X axis limits
        xlim(0,1000)
    
    # Print the plot in the pdf.
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
        proportion_cp2 = proportions[2],
        percentile_99 = cutoff_99,
        fc_99 = cutoff_99/means[2]
    )

}

# List of files to be processed.
filenames <- c("1A_EPCAM_measurements.txt",
               "2A_EPCAM_measurements.txt",
               "15A_EPCAM_measurements.txt",
               "20A_EPCAM_measurements.txt")

results <- filenames %>%
    # Map loops through all files and performs the tasks in each of them
    map_dfr(~ {  
        data <- read_tsv(.x) %>%  # load the data
            mutate(Image=as.character(str_replace_all(Image, "-|\\.tif$", "")))  # remove the extension of the file name in the Image column.

        # Create a PDF for the current file, this will create one pdf file with all the plots for each measurement file.
        pdf_filename <- str_replace(.x, "_measurements.txt", "_histograms.pdf")
        pdf(pdf_filename,  width = 390 / 25.4, height = 240 / 25.4)

        # Process each Image within the file and plot
        image_results <- data %>%
            group_by(Image) %>%
            do(plot_gmm_histogram(., unique(.$Image), .x))
        
        # Close the pdf only after all images were processed in the file
        dev.off()
        
        # Save the summary table that returned from the function.
        image_results
    })

# Save the summary table as a csv file.
results %>%
    write_tsv("gmm_table.txt")
