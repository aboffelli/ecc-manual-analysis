suppressMessages(suppressWarnings({
  source("renv/activate.R")

  input_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]

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

    df <- data.frame(
      pathologist = pathologist, tma = tma, sample = sample, paccs = paccs
    )
  }

  tidy <- input_file |>
    readLines() |>
    lapply(parse_line) |>
    Reduce(rbind, x = _)

  write.csv(tidy, output_file)
}))
