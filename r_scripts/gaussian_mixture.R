library(tidyverse)
library(gridExtra)
library(grid)
library(mclust)

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
        geom_histogram(aes(y = after_stat(density)), bins = 300, fill = "#619CFF", 
                       color = "#08519C", alpha=0.5) +
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

filenames <- c("2A_dab_measurements.txt",
               # "1A_measurements.txt",
               # "2A_measurements.txt",
               # "15A_measurements.txt",
               # "20A_measurements.txt",
               "2A_he_measurements.txt")

# dab_model <- read_tsv("2A_dab_model.tsv") %>% 
#     mutate(Source="DAB")
# merged <- bind_rows(he_model,dab_model) %>% 
#     mutate(Source=as.character(Source))

results <- filenames %>%
    map_dfr(~ {
        data <- read_tsv(.x) %>% 
            mutate(Image=as.character(str_remove(Image, "\\.tif$")))
        
        # Create a PDF for the current file
        pdf_filename <- str_replace(.x, "_measurements.txt", "_histograms.pdf")
        pdf(pdf_filename,  width = 390 / 25.4, height = 240 / 25.4)
        
        # Process each Image within the file and plot
        image_results <- data %>%
            group_by(Image) %>%
            do(plot_gmm_histogram(., unique(.$Image), .x))
    #     summary_table <- data %>% 
    #         group_by(Image) %>%
    #         summarise(
    #             # Fit a Gaussian Mixture Model with 2 components
    #             model = list(Mclust(`Area µm^2`, G = 2)),
    #             
    #             # Extract means for each component
    #             mean_cp1 = model[[1]]$parameters$mean[1],
    #             mean_cp2 = model[[1]]$parameters$mean[2],
    #             
    #             # Extract variances (diagonal of covariance matrix)
    #             variation_cp1 = model[[1]]$parameters$variance$sigmasq[1],
    #             variation_cp2 = model[[1]]$parameters$variance$sigmasq[2],
    #             
    #             # Extract proportions of each component
    #             proportion_cp1 = model[[1]]$parameters$pro[1],
    #             proportion_cp2 = model[[1]]$parameters$pro[2]
    #         ) %>%
    #         select(-model) %>% 
    #         mutate(Parent = str_remove(.x, "_measurements.txt"), 
    #                .before = Image)
    # 
    # return(summary_table)
        dev.off()
        image_results
    })
View(results)

# Plot for only one 
he_hist <- he_model %>% 
    ggplot(aes(x=`Area µm^2`)) +
    geom_histogram(bins = 50) +
    labs(title = "HE model")+
    theme_bw(); he_hist

dab_hist <- dab_model %>% 
    ggplot(aes(x=`Area µm^2`)) +
    geom_histogram(bins = 50) +
    labs(title = "DAB model")+
    theme_bw(); dab_hist

grid.arrange(he_hist, dab_hist, ncol = 2)

he_area <- he_model$`Area µm^2`
dab_area <- dab_model$`Area µm^2`


# Fit a gaussian mixture in the data
he_gmm_model <- Mclust(he_area, G=2)
summary(he_gmm_model)

# # Plot classification results
# plot(he_gmm_model, what = "classification")
# 
# # Plot BIC for model selection
# plot(he_gmm_model, what = "BIC")


# Fit a GMM to the data (univariate example)
data <- he_area
gmm_model <- Mclust(data, G=2)
summary(gmm_model)
component_means <- gmm_model$parameters$mean

# Create a sequence of points across the data range to plot the GMM line
x_seq <- seq(min(data), max(data), length.out = 1000)

# Calculate the GMM density for each component
gmm_density <- rep(0, length(x_seq))

# Loop through each component of the GMM
for (k in 1:gmm_model$G) {
    gmm_density <- gmm_density + gmm_model$parameters$pro[k] * 
        dnorm(x_seq, mean = gmm_model$parameters$mean[k], 
              sd = sqrt(gmm_model$parameters$variance$sigmasq[k]))
}

# Prepare the data frame for plotting
density_data <- data.frame(x = x_seq, y = gmm_density)

# Plot the histogram and overlay the GMM line
p1 <- ggplot(data.frame(x = data), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 300, fill = "#619CFF", 
                   color = "#08519C", alpha=0.5) +
    geom_line(data = density_data, aes(x = x, y = y), color = "black",
            linewidth = 1.2) +
    geom_vline(xintercept = component_means, linetype = c("dashed", "dotdash"), ,
              color = c("#A50026", "#E08214"), size = 2) +
    ggtitle("") +
    labs(x="Nuclear Area", y="Density") +
    theme_classic(base_size = 40)+
    xlim(0,1000); p1

ggsave("~/PACC/TMA/Manuscript/gmm_fit_and_means_histogram.tif", p1, width = 390, height=240, units = "mm")

# Normal distribution plot
data_mean <- mean(data)
data_sd <- sd(data)

# Create a sequence of points across the data range
x_seq <- seq(min(data), max(data), length.out = 1000)

# Calculate the normal distribution density at these points
normal_density <- dnorm(x_seq, mean = data_mean, sd = data_sd)

# Create a data frame for the normal distribution
normal_data <- data.frame(x = x_seq, y = normal_density)

# Plot the histogram with the normal distribution curve
p2 <- ggplot(data.frame(x = data), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 300, fill = "#619CFF",
                   color = "black") +
    geom_line(data = normal_data, aes(x = x, y = y), color = "#F8766D", 
              linewidth = 1.2) +
    ggtitle("Histogram with Normal Distribution Fit") +
    labs(x="Area", y="Density", title = "Normal distribution") +
    theme_bw()

grid.arrange(p2,p1, ncol = 2, 
             top = textGrob("DAB model area dist", 
                            gp = gpar(fontsize = 16, fontface = "bold")))
