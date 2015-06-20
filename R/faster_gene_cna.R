library(data.table)
library(plyr)
library(iGC)

library(doMC)
doMC::registerDoMC(cores=4)

sample_desc_pth <- system.file("extdata", "sample_desc.csv", package = "iGC")
sample_desc <- create_sample_desc(sample_desc_pth)


gain_th <- 0.4
loss_th <- -0.4

data("hg19DBNM", package = "iGC", envir = environment())
modified_hg19 <- hg19DBNM[
  , list(Start=min(.SD$Start), Stop=max(.SD$Stop)),
  by='Gene.Symbol,Chromosome'
]
# setkeyv(modified_hg19, c("Gene.Symbol", "Chromosome", "Start", "Stop"))
setkeyv(modified_hg19, c("Chromosome", "Start", "Stop"))
all_genes <- unique(modified_hg19[, .(Gene.Symbol)])
setkeyv(all_genes, c("Gene.Symbol"))

cna <- iGC:::read_cna(sample_desc$CNA_filepath[[1]])

library(doMC)
doMC::registerDoMC(cores=4)

joint_gene_cna <- aaply(
  sample_desc$CNA_filepath,
  1,
  function(cna_filepath) {
    cna <- iGC:::read_cna(cna_filepath)
    setkeyv(cna, c("Chromosome", "Start", "End"))

    g <- as.list(modified_hg19["BRCA2"])
    cna[
      Chromosome == g$Chromosome &
        Start <= g$Stop &
        End >= g$Start,
      .(Segment_Mean)
      ]
    cna[, gain_loss:=0]
    cna[Segment_Mean > gain_th, gain_loss:= 1]
    cna[Segment_Mean < loss_th, gain_loss:= -1]
    cna <- cna[, .(Chromosome, Start, End, gain_loss)]

    gene_overlapped <- foverlaps(
      cna, modified_hg19,
      by.x = c("Chromosome", "Start", "End"),
      by.y = c("Chromosome", "Start", "Stop"),
      type="any", nomatch=0L, mult = 'all'
    )[, .(Chromosome, Gene.Symbol, gain_loss)]

    gene_collapsed <- gene_overlapped[
      ,{n_gain <- sum(gain_loss == 1)
      n_norm <- sum(gain_loss == 0)
      n_loss <- sum(gain_loss == -1)
      if (n_gain == 0 && n_loss == 0) 0
      else ifelse(n_gain > n_loss, 1, -1)},
      keyby="Gene.Symbol"
      ]
    gene_collapsed[all_genes, .(V1),nomatch=NA][[1]]
  },
  .progress = "time",
  .parallel = TRUE
)
