suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(gtools)))


plot_gmm_histogram <- function(data, image_name, filename) {
    filename = str_remove(filename, "_measurements.txt")
    cat(paste("Processed: ", filename, image_name, "\n"))

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

# Access the input files
filenames <- snakemake@input

# Access output
output_file <- snakemake@output[[1]]

results <- filenames %>%
    map_dfr(~ {
        data <- suppressMessages(suppressWarnings(read_tsv(.x))) %>% 
            mutate(Image=as.character(str_replace_all(Image, "-|\\.tif$", "")))

        # Create a PDF for the current file
        filename <- basename(.x)
        pdf_filename <- str_replace(filename, 
                                    "_measurements.txt", 
                                    "_histograms.pdf")

        pdf(paste0("results/plots/", pdf_filename),  
            width = 390 / 25.4, 
            height = 240 / 25.4)
        
        # Process each Image within the file and plot
        image_results <- data %>%
            group_by(Image) %>%
            do(plot_gmm_histogram(., unique(.$Image), filename))
        dev.off()
        
        # Save the summary table that returned from the function.
        image_results
    })

results %>%
    write_tsv(output_file)
