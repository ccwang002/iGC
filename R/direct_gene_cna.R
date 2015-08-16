#' Load the existed CNA gain/loss based on gene location.
#'
#' This function aims to complement \code{\link{create_gene_cna}}. Instead of
#' mapping CNA records onto genes by genome reference, it reads the existed
#' column containing the gene each CNA lies on. Two functions share the same
#' interface but they have different requirement for the \code{read_fun}
#' implementation.
#'
#' @section Custom reader function: Similar to that of \link{create_gene_cna},
#'   the reader function takes the filepath as the first argument. It will
#'   return a data.table with at least two columns: \code{GENE} and
#'   \code{Segment_Mean} of type \code{character} and \code{numeric}
#'   respectively.
#'
#' @inheritParams create_gene_cna
#'
#' @return data.table of CNA gain/loss on each gene region for all samples,
#'   whose rows represent regions of genes and columns are sample names. First
#'   column \code{GENE} contains the corresponding gene names.
#'
#' @seealso \code{\link{create_gene_cna}}
#' @examples
#'
#' require(data.table)
#'
#' ## Create a CNA dataset that has been already mapped onto gene regions
#'
#' cna_geo_list = list(
#'     sample_A = data.table(
#'         GENE = c("TP53", "BRCA2"),
#'         Segment_Mean = c(1.05, -2.03)
#'     ),
#'     sample_B = data.table(
#'         GENE = c("TP53", "BRCA2", "NDPH1"),
#'         Segment_Mean = c(0.38, -1.71, 2.6)
#'     )
#' )
#' sample_desc <- data.table(
#'     Sample = paste("sample", c("A", "B"), sep = "_")
#' )
#' sample_desc$CNA_filepath <- sample_desc$Sample
#'
#'
#' ## Example code for reading
#'
#' read_cna_geo <- function(pth) {
#'     # For demonstration, file reading silently redirects
#'     # to list lookup
#'     cna_geo_list[[pth]]
#' }
#' gene_cna <- direct_gene_cna(
#'     sample_desc,
#'     read_fun = read_cna_geo, progress = FALSE
#' )
#' gene_cna
#'
#' @export
#' @import data.table
direct_gene_cna <- function(
    sample_desc,
    gain_threshold = log2(2.5) - 1,
    loss_threshold = log2(1.5) - 1,
    read_fun = NULL,
    progress = TRUE, progress_width = 48,
    parallel = FALSE,
    ...
) {
    # make the progress bar width smaller
    if (progress) {
        progress_opt <- "time"
        old_width_option <- options(width = progress_width)
    } else {
        progress_opt <- "none"
    }

    cna_bygene_list <- plyr::alply(
        sample_desc$CNA_filepath, 1,
        process_cna_per_sample_direct,
        gain_th = gain_threshold,
        loss_th = loss_threshold,
        read_fun = read_fun,
        ...,
        .progress = progress_opt, .parallel = parallel
    )
    all_genes <- unique(rbindlist(plyr::llply(
        cna_bygene_list, function(cna_dt) cna_dt[, .(GENE)]
    )))
    joint_gene_cna <- as.data.table(plyr::llply(
        cna_bygene_list,
        function(cna_bygene) {
            cna_bygene[all_genes, .(gain_loss), nomatch = NA][[1]]
        },
        .parallel = parallel
    ))
    # restore old width option
    if (progress) options(old_width_option)

    # restore sample and gene names
    all_samples <- sample_desc$Sample
    setnames(joint_gene_cna, seq_along(all_samples), all_samples)
    joint_gene_cna[, GENE := all_genes[[1]]]
    setcolorder(joint_gene_cna, c("GENE", all_samples))
    return(joint_gene_cna)
}

#' @import data.table
process_cna_per_sample_direct <- function(
    cna_filepath, gain_th, loss_th,
    read_fun = NULL,
    ...
) {
    if (is.null(read_fun)) read_fun <- read_cna_geo
    cna <- read_fun(cna_filepath, ...)
    setkeyv(cna, c("GENE"))
    # average multiple record for same gene
    cna_bygene <- cna[, list(cna_val = mean(Segment_Mean)), by = "GENE"]
    setkeyv(cna_bygene, c("GENE"))
    cna_bygene[, gain_loss := 0]
    cna_bygene[cna_val > gain_th, gain_loss := 1]
    cna_bygene[cna_val < loss_th, gain_loss := -1]
    cna_bygene
}

#' @import data.table
read_cna_geo <- function(cna_filepath) {
    dt <- fread(
        cna_filepath, sep = ',',
        colClasses = c("character", "character", "integer", "integer", "numeric")
    )
    setnames(dt, c("GENE", "Chromosome", "Start", "End", "Segment_Mean"))
    dt[, .(GENE, Segment_Mean)]
}
