#' @import data.table
#' @import plyr
#' @export
create_gene_exp <- function(sample_desc, read_fun = NULL, progress = TRUE, progress_width = 48, ...) {
  if (is.null(read_fun)) {
    read_fun <- read_gene_exp
  }
  ge_filepaths <- sample_desc$GE_filepath

  # make the progress bar width smaller
  if (progress) {
    progress_opt <- "time"
    old_width_option <- options(width = progress_width)
  } else {
    progress_opt <- "none"
  }
  raw_dfs <- alply(ge_filepaths, 1, read_fun, ..., .expand = FALSE, .progress = progress_opt)
  # restore old width option
  if (progress) options(old_width_option)

  df <- as.data.table(llply(raw_dfs, function(df){df[[2]]}))
  setnames(df, seq_len(ncol(df)), sample_desc$Sample)
  df$GENE <- raw_dfs[[1]][[1]]  # gene names from the parsed output of the first sample
  setcolorder(df, c("GENE", sample_desc$Sample))
  df
}

read_gene_exp <- function(ge_filepath, ...) {
  gene_exp <- read.table(
    ge_filepath,
    header = FALSE,
    skip = 2,
    na.strings = "null",
    colClasses = c("character", "double"),
    ...
  )
  # We have to read in as character and parse it as double. See bug report of
  # data.table on https://github.com/Rdatatable/data.table/issues/504
  # gene_exp <- fread(
  #   ge_filepath, header = FALSE, skip = 2,
  #   colClasses = c("character", "double"), na.strings=c("null")
  # )
  as.data.table(gene_exp)
}




