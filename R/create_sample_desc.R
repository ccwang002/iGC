#' Create sample description table containing all required inputs
#'
#' Each sample will have a unique name along with a pair of CNA and gene
#' expression file. This function generates a fixed format of sample description
#' by either reading an external CSV file or specifying them through separate
#' arugments with same order.
#'
#' @param sample_desc_filepath external sample description CSV file having at
#'   least these three columns: \code{Sample}, \code{CNA_filepath}, and
#'   \code{GE_filepath}. Note that the column names must be \emph{as is}.
#'
#' @param sample_names character vector of distinct sample names. Sample will be
#'   referenced by the given name through out the analysis process. They are
#'   expected to be valid names for R data.table's column names.
#'
#' @param cna_filepaths character vector of filepaths to CNA data.
#' @param ge_filepaths character vector of filepaths to gene expression data.
#' @param sample_root path to the root of sample data. If given, this path will
#'   be appended before all given filepaths.
#'
#' @return data.table of sample description having columns in order:
#'   \code{Sample}, \code{CNA_filepath}, and \code{GE_filepath}. Each row
#'   represents a sample by its unique sample name and the corresponding
#'   filepaths to CNA and gene expression data.
#'
#' @examples
#' ## Custom sample desc. by specifying separate arguments
#'
#' sample_names <- letters[1:5]
#' sample_desc <- create_sample_desc(
#'   sample_names = sample_names,
#'   cna_filepaths = file.path('cna', paste0(sample_names, '.csv')),
#'   ge_filepaths = file.path('ge', paste0(sample_names, '.txt'))
#' )
#' sample_desc
#'
#'
#' ## Prepend the file path with a root directory /path/to/sample
#'
#' create_sample_desc(
#'   sample_names = sample_desc$Sample,
#'   cna_filepaths = sample_desc$CNA_filepath,
#'   ge_filepaths = sample_desc$GE_filepath,
#'   sample_root = '/path/to/sample'
#' )
#'
#'
#' ## Create by reading the sample desc. CSV file
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
#' @note It's a good practice to specify absolute file paths. One could convert
#'   the relative file paths by passing \code{sample_root}.
#'
#'   If for special reason one does not want to specify one of the CNA of gene
#'   expression file paths, pass it with empty character vector of correct
#'   length, such as \code{rep('', num_samples)}.
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
    if (any(sapply(list(sample_names, cna_filepaths, ge_filepaths), is.null))) {
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
      sample_desc$CNA_filepath <- file.path(sample_root, sample_desc$CNA_filepath)
      sample_desc$GE_filepath <- file.path(sample_root, sample_desc$GE_filepath)
    }
  }
  sample_desc
}
