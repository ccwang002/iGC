#' Load and map CNA gain/loss by human gene location
#'
#' @param sample_desc data.table object created by \link{create_sample_desc}.
#' @param read_fun Custom reader function, see its own section for more detail.
#'
#' @import data.table
#' @export
create_gene_cna <- function(
  sample_desc,
  cna.gain.threshold = log(2.5, 2) - 1,
  cna.loss.threshold = log(1.5, 2) - 1,
  read_fun = NULL,
  progress = TRUE, progress_width = 48,
  ...
) {
  # TODO: don't ship our own hg19 DB
  data("hg19DBNM", package = "iGC", envir = environment())
  hg19DBNM <- data.table(hg19DBNM)

  # TODO: force no custom column names
  gene_wise_CNA <- unique(hg19DBNM[, .(Gene.Symbol)])
  orig_gene_order <- copy(gene_wise_CNA)
  setkey(gene_wise_CNA, Gene.Symbol) # add index

  # make the progress bar width smaller
  if (progress) {
    progress_opt <- "time"
    old_width_option <- options(width = progress_width)
  } else {
    progress_opt <- "none"
  }
  plyr::m_ply(
    sample_desc[, .(Sample, CNA_filepath)],
    read_cna_sample,
    gene_wise_CNA = gene_wise_CNA,
    gain_th = cna.gain.threshold,
    loss_th = cna.loss.threshold,
    progress = progress,
    read_fun = read_fun,
    ...,
    .progress = progress_opt
  )
  # restore old width option
  if (progress) options(old_width_option)

  setkey(gene_wise_CNA, NULL)  # remove the index
  # restore the order
  gene_wise_CNA[match(orig_gene_order[[1]], gene_wise_CNA[[1]])]
}

read_cna <- function(cna_filepath) {
  cna <- fread(cna_filepath, sep = '\t', header = TRUE)
  cna
}

read_cna_sample <- function(
  Sample, CNA_filepath,
  read_fun = NULL, gene_wise_CNA,
  gain_th, loss_th,
  progress = TRUE, ...
) {
  if (progress) cat("Sample:", Sample)
  if (is.null(read_fun)) {
    read_fun <- read_cna
  }
  cna <- read_fun(CNA_filepath, ...)
  gene_wise_CNA[, c(Sample):=0]
  for(cna_ix in seq_len(nrow(cna))) {
    cna_val <- cna[cna_ix, c(Segment_Mean)][[1]]
    if (cna_val > gain_th) {
      gol <- 1
    } else if(cna_val < loss_th) {
      gol <- -1
    } else {
      gol <- 0
      next   # early fail, don't care about no gain-or-loss records
    }
    cur_chr <- cna[cna_ix, .(Chromosome)][[1]]

    genes_on_same_chromosome <- hg19DBNM[Chromosome == cur_chr]
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

    # TODO: assure CNA recods cna[Start < End]
    genes_overlapped <- unique(genes_on_same_chromosome[
        Stop > cna[cna_ix, .(Start)][[1]] &
        Start < cna[cna_ix, .(End)][[1]],
      .(Gene.Symbol)][[1]])
    # TODO: judge by transcript, count 1 and -1, the more wins
    # Overwrite previous record if any
    gene_wise_CNA[genes_overlapped, c(Sample):=gol]
  }
}
