#' Create sample description table containing all required inputs
#'
#' Each sample will have a unique name along with a pair of CNA and gene
#' expression file. This function generates a table of sample descriptions by
#' either reading an external CSV file or specifying them through separate
#' arugments in same order.
#'
#' @param sample_desc_filepath external sample description CSV file having at
#'   least these three columns: \code{Sample}, \code{CNA_filepath}, and
#'   \code{GE_filepath}. Note that the column names must be given \emph{as is}.
#'
#' @param sample_names character vector of distinct sample names. Samples will
#'   be referenced by the given name through out the analysis process. They
#'   should be valid R data.table column names.
#'
#' @param cna_filepaths character vector of filepaths to CNA data.
#' @param ge_filepaths character vector of filepaths to gene expression data.
#' @param sample_root path to the root of sample data. If given, this path will
#'   be appended before all given filepaths.
#'
#' @return data.table of sample description having the following columns in
#'   order: \code{Sample}, \code{CNA_filepath}, and \code{GE_filepath}. Each row
#'   contains a sample's unique name and the corresponding filepaths to CNA and
#'   gene expression data.
#'
#' @examples
#' ## Custom sample description by specifying separate arguments
#'
#' sample_names <- letters[1:5]
#' sample_desc <- create_sample_desc(
#'     sample_names = sample_names,
#'     cna_filepaths = file.path('cna', paste0(sample_names, '.csv')),
#'     ge_filepaths = file.path('ge', paste0(sample_names, '.txt'))
#' )
#' sample_desc
#'
#'
#' ## Prepend the file path with a root directory /path/to/sample
#'
#' create_sample_desc(
#'     sample_names = sample_desc$Sample,
#'     cna_filepaths = sample_desc$CNA_filepath,
#'     ge_filepaths = sample_desc$GE_filepath,
#'     sample_root = '/path/to/sample'
#' )
#'
#'
#' ## Create by reading a sample description CSV file
#'
#' sample_desc_pth <- system.file("extdata", "sample_desc.csv", package = "iGC")
#' sample_desc <- create_sample_desc(sample_desc_pth)
#'
#'
#' \dontrun{
#' ## Read a external description and append the given file paths
#' create_sample_desc('/path/to/desc.csv', sample_root='/path/to/sample/root')
#' }
#'
#' @note One could convert the relative file paths into absolute paths by
#'   passing the root folder path to \code{sample_root}.
#'
#'   If for some special reasons, for example gene expression of all samples
#'   have been collected or the CNA records for each gene exist, but do not have
#'   the file paths to either CNA or gene expression data, pass it with empty
#'   character vector of correct length, such as \code{rep('', num_samples)}.
#'
#' @import data.table
#' @export
create_sample_desc <- function(
    sample_desc_filepath = NULL,
    sample_names = NULL,
    cna_filepaths = NULL,
    ge_filepaths = NULL,
    sample_root = NULL
) {
    if(is.null(sample_desc_filepath)) {
        if (any(sapply(list(
            sample_names, cna_filepaths, ge_filepaths
        ), is.null))) {
            stop("Unsufficient parameters (name, cna_pths, ge_pths) to create sample desc")
        }
        if (is.null(sample_root)) {
            sample_desc <- data.table(
                Sample = sample_names,
                CNA_filepath = cna_filepaths,
                GE_filepath = ge_filepaths
            )
        } else {
            sample_desc <- data.table(
                Sample = sample_names,
                CNA_filepath = file.path(sample_root, cna_filepaths),
                GE_filepath = file.path(sample_root, ge_filepaths)
            )
        }
    } else {
        sample_desc <- fread(sample_desc_filepath, sep = ",", header = TRUE)
        sample_root <- dirname(sample_desc_filepath)
        if (!is.null(sample_root)) {
            sample_desc$CNA_filepath <-
                file.path(sample_root, sample_desc$CNA_filepath)
            sample_desc$GE_filepath <-
                file.path(sample_root, sample_desc$GE_filepath)
        }
    }
    sample_desc
}
