source("renv/activate.R")
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
#suppressMessages(suppressWarnings(library(chisq.posthoc.test)))
options(scipen = 99)
options(digits = 22)

calculate_max <- function(x) {
    if (all(is.na(x))) {
        return(NA)
    } else {
        return(max(x, na.rm = TRUE))
    }
}

# Access the input files
patient_cores <- snakemake@input[["patient_cores"]]
patient_unique <- snakemake@input[["patient_ecc"]]
hif1a_file <- snakemake@input[["hif1a_data"]]

# Access the output files
hif1a_plot <- "results/plots/hif1a_activation_barplot.tif"
patient_hif <- snakemake@output[[1]]

data <- suppressMessages(suppressWarnings(read_tsv(patient_cores, na = "#NULL!"))) %>%
    mutate(hist_gra=as.integer(hist_gra)) %>% 
    filter(exkl == 0) %>% 
    droplevels() %>% 
    mutate(PARENT=paste0(PARENT,"A"))

hif1A_data <- suppressMessages(suppressWarnings(read_tsv(hif1a_file))) 

patient_data_with_hif <- merge(data, hif1A_data, 
                               by = c("ID", "PARENT"), all.x = TRUE) %>% 
    rename(HIF1A_NUCLEI=`FERNO_SWEBCG_HIF1A Nuclei`,
           HIF1A_NUCLEI_INT= `FERNO_SWEBCG_HIF1A INT Nuclei`)


PACC_HIF <- patient_data_with_hif %>% 
    filter(tumor==1) %>% 
    select(patnr, ID, PARENT, contains("HIF1A_")) %>% 
    group_by(patnr) %>% 
    summarize(HIF1A_NUCLEI = calculate_max(HIF1A_NUCLEI),
              HIF1A_NUCLEI_INT = calculate_max(HIF1A_NUCLEI_INT)) %>% 
    mutate(
        HIF1A_ACTIVATION = case_when(
        HIF1A_NUCLEI <= 1 & HIF1A_NUCLEI_INT <=1 ~ 0,
        HIF1A_NUCLEI == 1 & HIF1A_NUCLEI_INT >=2 ~ 1,
        HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT == 1 ~ 1,
        HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT >= 2 ~ 2)
        )

# Save a file with the HIF data.
patient_data <- suppressMessages(suppressWarnings(read_tsv(patient_unique,
                 na = "#NULL!"))) %>%
    filter(exkl == 0, tumor==1) %>%
    droplevels()

patient_data_with_hif_merged <- patient_data %>%
    mutate(PACC=as.integer(PACC),
           PACC_YN = if_else(PACC != 2, 0, 2)) %>%
    left_join(PACC_HIF, by="patnr") %>%
    relocate(PACC_YN, .after= PACC)

write.table(patient_data_with_hif_merged, 
            patient_hif, 
            sep = "\t", row.names = FALSE, quote = FALSE)

#Box plot with the HIF score
p1 <- patient_data_with_hif_merged %>% 
    group_by(PACC, HIF1A_ACTIVATION) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count),               # Total count within each PACC group
           percentage = (count / total_count) * 100) %>%
    drop_na() %>% 
    ggplot(aes(x=as.factor(PACC), y=percentage, 
               fill=as.factor(HIF1A_ACTIVATION)))+
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "",
         x = "ECC Score", 
         y = "Percentage (%)", 
         fill = "HIF1A Activation\nScore") +
    scale_fill_brewer(palette = "Dark2")+
    theme_classic(base_size = 20) +
    theme(legend.position = c(0.80, 0.90),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18))


ggsave(hif1a_plot, p1, width=150, height=150, units = 'mm')




# contingency_table <- table(PACC_HIF$PACC_YN, PACC_HIF$HIF1A_ACTIVATION)
# print("Contingency Table:")
# print(contingency_table)      

# # Step 3: Calculate the expected matrix
# expected_matrix <- chisq.test(contingency_table)$expected
# print("Expected Matrix:")
# print(expected_matrix)

# # Step 4: Perform chi-square test
# result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
# result <- chisq.test(contingency_table)
# print("Chi-square Test Result:")
# print(result)

# posthoc_result <- chisq.posthoc.test(contingency_table, method = "bonferroni") %>% 
#     rename(HIF1A_Negative="0", HIF1A_Low="1", HIF1A_High="2") %>% 
#     mutate(Dimension=paste("PACC", Dimension))
# write_tsv(file="~/PACC/TMA/Results/HIF1A/post_hoc_pacc.tab", x=posthoc_result)

# ## Check HIF1A with hist_gra ----
# HIST_HIF <- patient_data_with_hif %>% 
#     filter(tumor==1) %>% 
#     select(patnr, hist_gra, contains("HIF1A_")) %>% 
#     group_by(patnr) %>% 
#     #mutate(hist_gra= if_else(PACC == "NA", NA, PACC)) %>%  
#     summarize(hist_gra=calculate_max(hist_gra),
#               HIF1A_NUCLEI = calculate_max(HIF1A_NUCLEI),
#               HIF1A_NUCLEI_INT = calculate_max(HIF1A_NUCLEI_INT)) %>% 
#     mutate(HIF1A_ACTIVATION = case_when(
#         HIF1A_NUCLEI <= 1 & HIF1A_NUCLEI_INT <=1 ~ 0,
#         HIF1A_NUCLEI == 1 & HIF1A_NUCLEI_INT >=2 ~ 1,
#         HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT == 1 ~ 1,
#         HIF1A_NUCLEI >= 2 & HIF1A_NUCLEI_INT >= 2 ~ 2),
#         HIST_GRA_YN= if_else(hist_gra != 3, 1, 3))


# contingency_table <- table(HIST_HIF$HIST_GRA_YN, HIST_HIF$HIF1A_ACTIVATION)
# print("Contingency Table:")
# print(contingency_table)      

# # Step 3: Calculate the expected matrix
# expected_matrix <- chisq.test(contingency_table)$expected
# print("Expected Matrix:")
# print(expected_matrix)

# # Step 4: Perform chi-square test
# result <- chisq.test(contingency_table)
# print("Chi-square Test Result:")
# print(result)

# posthoc_result <- chisq.posthoc.test(contingency_table, method = "bonferroni") %>% 
#     rename(HIF1A_Negative="0", HIF1A_high="2") %>% 
#     mutate(Dimension=paste("Grade", Dimension))
# write_tsv(file="~/PACC/TMA/Results/HIF1A/post_hoc_hist.tab", x=posthoc_result)


# ## PACC and Histology grade
# PACC_HIST <- patient_data_with_hif %>% 
#     filter(tumor==1) %>% 
#     select(patnr, PACC, hist_gra) %>% 
#     group_by(patnr) %>% 
#     mutate(PACC= if_else(PACC == "NA", NA, PACC)) %>%  
#     summarize(PACC=calculate_max(PACC),
#               hist_gra=calculate_max(hist_gra)) %>% 
#     mutate(HIST_GRA_YN= if_else(hist_gra != 3, 1, 3),
#            PACC_YN = if_else(PACC != 2, 0, 2))

# contingency_table <- table(PACC_HIST$HIST_GRA_YN, PACC_HIST$PACC_YN)
# print(contingency_table)  

# expected_matrix <- chisq.test(contingency_table)$expected
# print("Expected Matrix:")
# print(expected_matrix)

# # Step 4: Perform chi-square test
# result <- chisq.test(contingency_table)
# print("Chi-square Test Result:")
# print(result)

# posthoc_result <- chisq.posthoc.test(contingency_table, method = "bonferroni") %>% 
#     rename(PACC_0="0", PACC_2="2") %>% 
#     mutate(Dimension=paste("Grade", Dimension))


# print(patient_data_with_hif %>%
#         filter(PACC == 2, HIF1A_NUCLEI >= 2, HIF1A_NUCLEI_INT >= 2) %>%
#         select(PARENT, ID) %>%
#         mutate(
#             PARENT = format(PARENT, justify = "left"),
#             ID = format(ID, justify = "left")
#         ), row.names = F)
