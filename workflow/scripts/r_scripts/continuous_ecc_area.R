if (!interactive()) {
  stardist_measurements <- snakemake@input$stardist_rois
  pathologist_measurements <- snakemake@input$merged_ecc_sizes
  tma_map <- snakemake@input$tma_map
  output_file <- snakemake@output[[1]]
}

if (interactive()) {
  stardist_measurements <- paste0("./resources/QuPathDetection/", c("1", "2", "15", "20"), "A_EPCAM_measurements.txt")
  pathologist_measurements <- "./results/merged_ecc_sizes.txt"
  tma_map <- "./resources/Koppling_PathXL_till_TMA-nr.csv"
}

# Step 1: Get a range of nuclei area that pathologists have noted as 'big'
#
# ==============================================================================
#
# Start by parsing the file
#
# Lines are tab delimited, but the last column contains (possibly) several
# values (also tab delimited), so the file is 'ragged'.
#
# The first three columns are static - the first is pathologist (either
# Pathologist 1 or Pathologist 2). The second is the name of the TMA. The third
# is the coordinate of the sample in the TMA.
#
# We capture the first three and then deal with the rest as a vector. The rest
# of the columns will be recycled to the length of this vector
#
# If there are no values, return a data.frame that has NA for the paccs col
parse_line <- function(line) {
  values <- strsplit(line, "\t")[[1]]
  pathologist <- values[1]
  tma <- values[2]
  sample <- values[3]
  paccs <- values[-(1:3)]
  paccs <- paccs[paccs != ""] # removes trailing \t\t\t

  if (length(paccs) == 0) {
    paccs <- NA_real_
  } else {
    # ECC areas are reported with decimal comma - must replace before coercion
    paccs <- sub(",", ".", paccs, fixed = TRUE) |>
      as.numeric()
  }

  data.frame(
    pathologist = pathologist, tma = tma, sample = sample, paccs = paccs
  )
}

pathologist_tidy <- pathologist_measurements |>
  readLines() |>
  lapply(parse_line) |>
  Reduce(rbind, x = _)

# Step 2: Use this range to filter the ROIs provided by StarDist
# ==============================================================================

# Start by combining the stardist files
combine_stardist_files <- function(paths) {
  slide_names <- sapply(
    regmatches(paths, gregexpr("QuPathDetection/\\d+", paths)),
    \(x) substr(x, nchar("QuPathDetection/") + 1, nchar(x))
  )
  df_list <- lapply(paths, read.delim)

  # Slide is only encoded in file name - add explicitly as column
  mapply(
    \(x, y) {
      x$slide <- y
      x
    }, df_list, slide_names, SIMPLIFY = FALSE
  ) |>
    Reduce(rbind, x = _)
}

stardist_df <- combine_stardist_files(stardist_measurements)

# Filter out cells based on range given by pathologist-selected big cells
pacc_area_range <- range(pathologist_tidy$paccs, na.rm = TRUE)
big_cells <- stardist_df[stardist_df$`Area.µm.2` > pacc_area_range[1], ]
big_cells <- big_cells[big_cells$`Area.µm.2` < pacc_area_range[2], ]

# Step 3: Get a per-sample metric
#
# ==============================================================================
#
# Metric is based on the area above the range lower end. This way the metric
# will start at 0, which hopefully makes sample that have just a few PACCs
# comparable to those that have none
big_cells$area_above_min <- big_cells$`Area.µm.2` - pacc_area_range[1]

# Calculate sum of the metric per sample. This is to mimic the idea of using the
# number of PACCs, but it's weighted based on area now.
big_cells$group <- paste0(big_cells$Image, "_", big_cells$slide)
metric_list <- tapply(big_cells, ~ group, \(x) sum(x$area_above_min, na.rm = TRUE), simplify = FALSE)

metric_df <- metric_list |>
  names() |>
  strsplit("_") |>
  lapply(\(x) data.frame(Image = x[1], slide = x[2])) |>
  Reduce(rbind, x = _)

metric_df$metric <- unlist(metric_list)

# Step 4: Get a common ID for joining
#
# ==============================================================================
#
# I need to be able to join these data with the rest of the pipeline. One common
# column used is the TMAnr column, which we can get from the pathology data

tma_map <- read.delim(tma_map)

# Modify metric_df cols to match the format in tma_map
metric_df$Image <- regmatches(metric_df$Image, gregexpr("^[^\\.]*", metric_df$Image)) |>
  unlist() |>
  gsub("-", "", x = _)

metric_df$slide <- paste("Glas", metric_df$slide)

final <- merge(
  tma_map, metric_df,
  by.x = c("ID", "PARENT"), by.y = c("Image", "slide"),
  all.x = TRUE
)

write.csv(final, output_file)
