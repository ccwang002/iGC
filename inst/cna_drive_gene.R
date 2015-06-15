find_cna_driven_gene <- function(
  gene_cna, gene_exp,
  gain_ratio, loss_ratio
) {


}

library(data.table)
library(plyr)

all_samples <- colnames(gene_cna)[-1]

# get shared genes
shared_genes <- intersect(
  gene_cna[, .(Gene.Symbol)][[1]],
  gene_exp[, .(GENE)][[1]]
)


# transpose so columns to be gene-wise, easier to compute
setkeyv(gene_cna, c("Gene.Symbol"))
gene_cna_t <- t(gene_cna[shared_genes, all_samples, with=FALSE])
# see data.table's setattr doc page.
# The following equals to
#     dinames(gene_cna_t) <- list(NULL, shared_genes)
# setattr prevent from creating unneccessary object copies.
setattr(gene_cna_t, 'dimnames', list(NULL, shared_genes))

gol_ratio_table <- as.data.table(aaply(
  gene_cna_t,
  2,
  function(one_gene) {
    c(Gain = sum(one_gene == 1),
      Loss = sum(one_gene == -1),
      Normal = sum(one_gene == 0))},
  .progress = 'time'
) / length(all_samples))

setkeyv(gene_exp, c("GENE"))
gene_exp_t <- t(gene_exp[shared_genes, all_samples, with=FALSE])
setattr(gene_exp_t, 'dimnames', list(NULL, shared_genes))
