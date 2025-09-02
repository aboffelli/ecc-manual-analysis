source("renv/activate.R")
library(readr)

if (!interactive()) {
  all_data_path <- snakemake@input$all_data
  ecc_metrics_path <- snakemake@input$ecc_metric_per_core
  output_file_path <- snakemake@output[[1]]
}

# Mainly for development
if (interactive()) {
  all_data_path <- "./results/all_patient_unique_with_ecc_hif1a_gex.txt"
  ecc_metrics_path <- "./results/big_cell_rois.txt"
  output_file_path <- "./results/all_data_with_continuous_ecc.txt"
}

# Read in ecc_metrics
ecc_metrics <- read.csv(ecc_metrics_path)

# TMAnr for joining, metric for downstream calculation - rest is unneeded
ecc_metrics <- ecc_metrics[c("TMAnr", "metric")]

# We're using readr::read_tsv here instead of read.delim because *something* is
# going on in this file where a print to console causes emacs to crash after
# printing thousands and thousands of lines.
#
# Upon read-in with read.delim, there's a complaint about an EOF within a string
# or something - that might be the issue.
#
# Read in all_data
all_data <- read_tsv(all_data_path)

# Join ecc_metrics and all_data
both <- merge(
  all_data, ecc_metrics, by.x = "tma_nr", by.y = "TMAnr", all.x = TRUE
)

write.csv(both, output_file_path)
