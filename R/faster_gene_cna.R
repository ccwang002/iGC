#' Load and map CNA gain/loss by human gene location
#'
#' @param sample_desc data.table object created by \link{create_sample_desc}.
#' @param read_fun Custom reader function, see its own section for more detail.
#'
#' @import data.table
#' @export
faster_gene_cna <- function(
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

  joint_gene_cna <- as.data.table(alply(
    sample_desc$CNA_filepath, 1,
    process_cna_per_sample,
    gain_th = gain_threshold,
    loss_th = loss_threshold,
    read_fun = read_fun,
    gene_db = modified_hg19,
    all_genes = all_genes,
    ...,
    .progress = "time", .parallel = parallel
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


read_cna_tcga <- function(cna_filepath) {
  cna <- fread(cna_filepath, sep = '\t', header = TRUE)
  cna
}



