#' Load and map CNA gain/loss onto human gene location by genome reference
#'
#' The function reads through in all sample CNA data given by the sample
#' description \code{sample_desc} and returns a joint CNA gain/loss table based
#' on gene regions across samples.
#'
#' A gene is considered to have CNA gain if the overlapped CNA record expression
#' is higher than the given threshold. Similarly, a gene is considered CNA loss
#' if the overlapped CNA record is lower than the given threshold. If multiple
#' CNA records map onto the same gene region with both gain and loss, the
#' majority wins. If none of the records map to the gene, NA is given.
#'
#' By default it assumes the data to be of TCGA level 3 file format. For other
#' data formats (e.g. raw data or other experiments from GEO), one should
#' implement a custom reader function that accepts the filepath as the first
#' argument. See section \emph{Custom reader function} for full specification.
#'
#' Currently the package ships a custom genome reference hg19, \link{hg19DBNM},
#' for gene region look up. Each gene's region is defined by the widest splicing
#' form it has in NCBI curated records. The defined region includes intron
#' regions. This limitation may change in the future.
#'
#' @section Custom reader function: Custom reader function is given by
#'   \code{read_fun = your_reader_fun}. It takes the filepath to CNA data as the
#'   first argument and returns a data.table with at least the following four
#'   columns: \code{Chromosome}, \code{Start}, \code{End}, and
#'   \code{Segment_Mean} of type character, integer, integer and numeric
#'   respectively.
#'
#'   Rest arguments of \code{create_gene_cna(...)} will be passed to this reader
#'   function.
#'
#'   Note: all string-like columns should \strong{NOT} be of type \code{factor}.
#'   Remember to set \code{stringsAsFactors = FALSE}.
#'
#' @param sample_desc \link[data.table]{data.table} object created by
#'   \link{create_sample_desc}.
#' @param gain_threshold CNA expression above this will be considered as gain
#'   region. By default \eqn{\log_2{2.5} - 1}
#' @param loss_threshold CNA expression below this will be considered as loss
#'   region. By default \eqn{\log_2{1.5} - 1}
#' @param read_fun Custom reader function, see its own section for more detail.
#' @param progress Whether to display a progress bar. By default \code{TRUE}.
#' @param progress_width The text width of the shown progress bar. By default is
#'   48 chars wide.
#' @param parallel Enable parallelism by plyr. One has to specify a parallel
#'   engine beforehand. See example for more information.
#' @param ... Arguments passed to the custom reader function specified in
#'   \code{read_fun}.
#'
#' @return data.table of CNA gain/loss on each gene region for all samples,
#'   whose rows represent regions of genes and columns represent sample names.
#'   First column \code{GENE} contains the corresponding gene names.
#'
#' @seealso \code{\link[utils]{read.table}} and \code{\link[data.table]{fread}}
#'   for custom reader function implementation; \code{\link{create_sample_desc}}
#'   for creating sample description. If the gene information already exists in
#'   the data, try \link{direct_gene_cna} to skip the genome reference lookup.
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
#'
#' ## Define custom reader function for TCGA level 3 gene exp. data
#'
#' my_cna_reader <- function(cna_filepath) {
#'     cna <- data.table::fread(cna_filepath, sep = '\t', header = TRUE)
#'     data.table::setnames(
#'         cna,
#'         c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
#'     )
#'     # pick only the needed columns
#'     cna[, .(Chromosome, Start, End, Segment_Mean)]
#' }
#'
#'
#' ## Read all samples' CNA data and combined as a single table
#'
#' gene_cna <- create_gene_cna(
#'     sample_desc,
#'     gain_threshold = log2(2.3) - 1, loss_threshold = log2(1.7) - 1,
#'     read_fun = my_cna_reader,
#' )
#' gene_cna[GENE %in% c("BRCA2", "TP53", "SEMA5A"), ]
#'
#'
#' \dontrun{
#' ## To boost the speed, utilize parallelization
#'
#' doMC::registerDoMC(4)  # number of CPU cores
#' gene_cna <- create_gene_cna(
#'     sample_desc,
#'     gain_threshold = log2(2.3) - 1, loss_threshold = log2(1.7) - 1,
#'     read_fun = my_cna_reader,
#'     parallel = TRUE
#' )
#' }
#' @import data.table
#' @export
create_gene_cna <- function(
    sample_desc,
    gain_threshold = log2(2.5) - 1,
    loss_threshold = log2(1.5) - 1,
    read_fun = NULL,
    progress = TRUE, progress_width = 48,
    parallel = FALSE,
    ...
) {
    # TODO: don't ship our own hg19 DB, make other ref customizable
    data("hg19DBNM", package = "iGC", envir = environment())

    # select the largest range for each gene (on each chromosome)
    modified_hg19 <- hg19DBNM[
        , list(Start=min(.SD$Start), Stop=max(.SD$Stop)),
        by='Gene.Symbol,Chromosome'
        ]
    setkeyv(modified_hg19, c("Chromosome", "Start", "Stop"))

    # TODO: actually we can use other identifier to prevent collision
    # select unique gene names
    all_genes <- unique(modified_hg19[, .(Gene.Symbol)])
    setkeyv(all_genes, c("Gene.Symbol"))

    # make the progress bar width smaller
    if (progress) {
        progress_opt <- "time"
        old_width_option <- options(width = progress_width)
    } else {
        progress_opt <- "none"
    }

    joint_gene_cna <- as.data.table(plyr::alply(
        sample_desc$CNA_filepath, 1,
        process_cna_per_sample,
        gain_th = gain_threshold,
        loss_th = loss_threshold,
        read_fun = read_fun,
        gene_db = modified_hg19,
        all_genes = all_genes,
        ...,
        .progress = progress_opt, .parallel = parallel
    ))

    # restore old width option
    if (progress) options(old_width_option)

    # restore sample and gene names
    all_samples <- sample_desc$Sample
    setnames(joint_gene_cna, seq_along(all_samples), all_samples)
    joint_gene_cna[, GENE:=all_genes[[1]]]
    setcolorder(joint_gene_cna, c("GENE", all_samples))
    return(joint_gene_cna)
}


#' @import data.table
process_cna_per_sample <- function(
    cna_filepath, gain_th, loss_th,
    read_fun = NULL, gene_db = NULL, all_genes = NULL,
    ...
) {
    if (is.null(read_fun)) read_fun <- read_cna_tcga

    cna <- read_fun(cna_filepath, ...)
    setkeyv(cna, c("Chromosome", "Start", "End"))
    cna[, gain_loss:=0]
    cna[Segment_Mean > gain_th, gain_loss:= 1]
    cna[Segment_Mean < loss_th, gain_loss:= -1]
    cna <- cna[, .(Chromosome, Start, End, gain_loss)]

    # Overlap table, implemented by data.table::foverlaps
    # CNA    a        ++++++++++++++
    # Gene   A     >>>>>-->>>>>>>->>>>
    #   exon 1     T&T                                T
    #   exon 2            T&T                         T
    #   exon 3                    T&T                 T
    # Gene   B                       <<<<---<<<<<<<
    #   exon 1                       T&F              F
    #   exon 2                              T&F       F
    # Gene   C     >>>>>>>>>>>>>>>>>>>>>
    #   exon 1     T&T                                T
    # Gene   D <<<
    #   exon 1 F&T                                    F
    gene_overlapped <- foverlaps(
        cna, gene_db,
        by.x = c("Chromosome", "Start", "End"),
        by.y = c("Chromosome", "Start", "Stop"),
        type="any", nomatch=0L, mult = 'all'
    )[, .(Chromosome, Gene.Symbol, gain_loss)]

    gene_collapsed <- gene_overlapped[
        ,list(cur_sample = judge_gain_or_loss(gain_loss)),
        keyby="Gene.Symbol"
        ]
    return(gene_collapsed[all_genes, .(cur_sample), nomatch=NA][[1]])
}


judge_gain_or_loss <- function(gain_loss) {
    n_gain <- sum(gain_loss == 1)
    n_norm <- sum(gain_loss == 0)
    n_loss <- sum(gain_loss == -1)
    if (n_gain == 0 && n_loss == 0) {
        return(ifelse(n_norm > 0, 0, NA))
    } else {
        return(ifelse(n_gain > n_loss, 1, -1))
    }
}


#' @import data.table
read_cna_tcga <- function(cna_filepath) {
    cna <- fread(cna_filepath, sep = '\t', header = TRUE)
    cna
}
