#' @import data.table
#' @export
create_sample_desc <- function(
  sample_desc_filepath = NULL,
  sample_names = NULL,
  cna_filepaths = NULL,
  ge_filepaths = NULL,
  relative_path = TRUE
) {
  if(is.null(sample_desc_filepath)) {
    if (any(sapply(list(sample_names, cna_filepaths, ge_filepaths), is.null))) {
      stop("Unsufficient parameters (name, cna_pths, ge_pths) to create sample desc")
    }
    sample_desc <- data.table(
      Sample = sample_names,
      CNA_filepath = cna_filepaths,
      GE_filepath = ge_filepaths
    )
  } else {
    sample_desc <- fread(sample_desc_filepath, sep = ",", header = TRUE)
    sample_root <- dirname(sample_desc_filepath)
    if (relative_path) {
      sample_desc$CNA_filepath <- file.path(sample_root, sample_desc$CNA_filepath)
      sample_desc$GE_filepath <- file.path(sample_root, sample_desc$GE_filepath)
    }
  }
  sample_desc
}
