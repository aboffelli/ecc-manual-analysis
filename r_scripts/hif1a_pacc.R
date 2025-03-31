library(tidyverse)
library(RColorBrewer)
library(chisq.posthoc.test)
options(scipen = 99)
options(digits = 22)

# Function to transform the percentages labels in integers.
map_percent_to_numeric <- function(x) {
    case_when(
        x == "0%" ~ 0,
        x == "1-10%" ~ 1,
        x == "11-50%" ~ 2,
        x == ">50%" ~ 3,
        x == "No score" ~ NA_real_,
        TRUE ~ NA_real_
    )
}

# Function to retrieve the highest number
calculate_max <- function(x) {
    if (all(is.na(x))) {
        return(NA)
    } else {
        return(max(x, na.rm = TRUE))
    }
}

# Set the directory variable to avoid retyping it
dir <- "~/PACC/TMA/Results/"

# Load the patient data that contains the TMA numbers
data <- read_tsv(paste0(dir,
                        "all_patients_with_TMA_PACC_new.csv"),na = "#NULL!") %>%
    mutate(hist_gra=as.integer(hist_gra)) %>%  # Make the grade into integer.
    filter(exkl == 0) %>%  # Filter only patients that were used in the original study.
    droplevels() %>% 
    mutate(PARENT=paste0(PARENT,"A"))  # Add an "A" to the TMA number to match the other file.

# Load the hif scores file
hif1A_data <- read_tsv("~/PACC/TMA/Data/Tables/ImunoStainings/hif_scores_merged.csv") 

# Merge the data based on the Core ID and TMA number (PARENT).
patient_data_with_hif <- merge(data, hif1A_data, 
                               by = c("ID", "PARENT"), all.x = TRUE) %>% 

    # Rename columns to make our life easier.
    rename(HIF1A_NUCLEI=`FERNO_SWEBCG_HIF1A Nuclei`,
           HIF1A_NUCLEI_INT= `FERNO_SWEBCG_HIF1A INT Nuclei`)


PACC_HIF <- patient_data_with_hif %>% 
    
    # Filter only primary tumor.
    filter(tumor==1) %>%

    # Select only patient number, core ID, TMA number and the HIF scores. 
    select(patnr, ID, PARENT, contains("HIF1A_")) %>% 

    group_by(patnr) %>% 

    # Calculate the max score for each patient using the custom function.
    summarize(HIF1A_NUCLEI = calculate_max(HIF1A_NUCLEI),
              HIF1A_NUCLEI_INT = calculate_max(HIF1A_NUCLEI_INT)) %>% 

    # Merge the HIF score using area and intensity scores
    mutate(
        HIF1A_ACTIVATION = case_when(

        HIF1A_NUCLEI <= 1 & HIF1A_NUCLEI_INT <=1 ~ 0,   # area and intesity lower or equal to 1-10%  
        HIF1A_NUCLEI == 1 & HIF1A_NUCLEI_INT >=2 ~ 1,   # area between 1-10% and intensity equal or over than 11-50%
        HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT == 1 ~ 1,  # area equal or over than 11-50% and intesity between 1-10%
        HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT >= 2 ~ 2)  # area and intesity equal or over than 11-50%
        )

## Plot ----

# Bar plot to visualize the percentage of patients are in each group of HIF1A Activation
p1 <- PACC_HIF %>% 
    group_by(PACC_YN, HIF1A_ACTIVATION) %>%
    summarise(count = n()) %>%  # Count the number of patients
    mutate(total_count = sum(count),    # Total count within each PACC_YN group
           percentage = (count / total_count) * 100) %>%    #Calculate percentage
    drop_na() %>%   # Remove NAs

    # Barplot of the percentages
    ggplot(aes(x=as.factor(PACC_YN), y=percentage, 
               fill=as.factor(HIF1A_ACTIVATION)))+  # Fill by HIF1A group
    
    geom_bar(stat = "identity", position = "dodge") +   # Bars side by side. 
    
    # Change the labels
    labs(title = "",
         x = "ECC Score", 
         y = "Percentage (%)", 
         fill = "HIF1A Activation\nScore") +
        
    scale_fill_brewer(palette = "Dark2")+   # New pallete that is colorblind friendly
    
    theme_classic(base_size = 20) +     # Change the theme and increase font size
    
    # Put the legend inside the plot and increase the font.
    theme(legend.position = c(0.80, 0.90),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18))

# Save the tif image
ggsave("~/PACC/TMA/Manuscript/hif1a_activation_barplot.tif", p1, width=150, height=150, units = 'mm')

# Load the patient file with unique patients.
patient_data <- read_tsv(paste0(dir,
                        "all_patient_unique_with_pacc_full_score_new.csv"),
                 na = "#NULL!") %>%

    # exclude patients removed from the original study and keep only primary tumor
    filter(exkl == 0, tumor==1) %>% 
    droplevels()

# Merge the HIF1A data in the patient data
patient_data_with_hif_merged <- patient_data %>%
    mutate(PACC=as.integer(PACC),
           PACC_YN = if_else(PACC != 2, 0, 2)) %>%

    left_join(PACC_HIF, by="patnr") %>% # Use the patient number to merge

    relocate(PACC_YN, .after= PACC) # Move PACC_YN close to PACC

# Save it as a tab delimited file
write.table(patient_data_with_hif_merged,
            paste0(dir,
                   "all_patients_with_TMA_PACC_HIF.csv"),
            sep = "\t", row.names = FALSE, quote = FALSE)