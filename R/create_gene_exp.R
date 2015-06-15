#' @import data.table
#' @import plyr
#' @export
create_gene_exp <- function(sample_desc, read_fun = NULL, ...) {
  if (is.null(read_fun)) {
    read_fun <- read_gene_exp
  }
  # use first_ge to obtain gene exp list
  ge_filepaths <- sample_desc$GE_filepath
  first_ge <- read_fun(ge_filepaths[1])
  raw_dfs <- alply(ge_filepaths, 1, read_fun, ..., .expand = FALSE, .progress = "time")
  # TODO: remove magic [[2]]
  df <- as.data.table(llply(raw_dfs, function(df){df[[2]]}))
  setnames(df, seq_len(ncol(df)), sample_desc$Sample)
  df$GENE <- first_ge[[1]]
  # row.names(df) <- ge[[1]]
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
  gene_exp
}




