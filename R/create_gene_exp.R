#' Create an joint gene expression table of all samples
#'
#' The function reads in all gene expression data given by the sample
#' description \code{sample_desc} and return a joint expression table of all
#' samples.
#'
#' By default it assumes the data to be of TCGA level 3 file format. However,
#' nearly all real world data fail to have the same format as TCGA. In this
#' case, one needs to tell the function how to parse the data by implementing a
#' custom reader function that accepts the filepath as the first argument. See
#' Detail section for full specification. The function naively concatenates all
#' return expression \emph{as if all gene expressions are stated in the same
#' gene order} as columns in a new data.table.
#'
#' @section Custom reader function: Custom reader function is given by
#'   \code{read_fun = your_reader_fun}. It takes the filepath as the first
#'   argument and return a data.table with the first two columns being
#'   \code{GENE} and \code{Expression} of type character and double.
#'
#'   The output joint gene expression table has first column \code{GENE} store
#'   the gene name, which are are determined by the first sample being
#'   evaluated.
#'
#'   Rest arguments of \code{create_gene_exp(...)} will be passed to this reader
#'   function.
#'
#'   Note: all string-like columns should \strong{NOT} be of type \code{factor}.
#'   Remember to set \code{stringsAsFactors = FALSE}.
#'
#' @param sample_desc data.table object created by \link{create_sample_desc}.
#' @param read_fun Custom reader function, see its own section for more detail.
#' @param progress Whether to display a progress bar. By default \code{TRUE}.
#' @param progress_width The text width of the shown progress bar. By default is
#'   48 chars wide.
#' @param ... Arguments passed to the custom reader function specified in
#'   \code{read_fun}.
#'
#' @return data.table of all samples gene expression, whose rows are gene
#'   expression and columns are sample names. First column \code{GENE} contains
#'   the corresponding gene names.
#'
#' @seealso \code{\link[utils]{read.table}} and \code{\link[data.table]{fread}}
#'   for custom reader function implementation; \code{\link{create_sample_desc}}
#'   for creating sample description.
#'
#' @examples
#' ## Use first three samples of the builtin dataset
#'
#' sample_root <- system.file("extdata", package = "iGC")
#' sample_desc_pth <- file.path(sample_root, "sample_desc.csv")
#' sample_desc <- create_sample_desc(
#'     sample_desc_pth, sample_root=sample_root
#' )[1:3]
#'
#' ## Define custom reader function for TCGA level 3 data
#' my_gene_exp_reader <- function(ge_filepath) {
#'     gene_exp <- read.table(
#'         ge_filepath,
#'         header = FALSE, skip = 2,
#'         na.strings = "null",
#'         colClasses = c("character", "double")
#'     )
#'     dt <- data.table::as.data.table(gene_exp)
#'     data.table::setnames(dt, c("GENE", "Expression"))
#' }
#' gene_exp <- create_gene_exp(
#'     sample_desc,
#'     read_fun = my_gene_exp_reader,
#'     progress_width = 60
#' )
#' gene_exp[1:5]
#'
#' @note The function assumes row order for all samples' gene expressions are the
#'   same.
#' @import data.table
#' @import plyr
#' @export
create_gene_exp <- function(
    sample_desc, read_fun = NULL,
    progress = TRUE, progress_width = 48, ...
) {
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
    raw_dfs <- alply(
        ge_filepaths, 1, read_fun, ...,
        .expand = FALSE, .progress = progress_opt
    )
    # restore old width option
    if (progress) options(old_width_option)

    df <- as.data.table(llply(raw_dfs, function(df){df[[2]]}))
    setnames(df, seq_len(ncol(df)), sample_desc$Sample)
    # gene names from the parsed output of the first sample
    df$GENE <- raw_dfs[[1]][[1]]
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
