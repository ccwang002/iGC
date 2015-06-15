#' @import data.table
#' @import plyr
#' @export
find_cna_driven_gene <- function(
  gene_cna, gene_exp,
  gain_ratio = 0.2, loss_ratio = 0.2
) {
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
  gol_ratio_table[, GENE:=shared_genes]
  setkeyv(gol_ratio_table, "GENE")

  setkeyv(gene_exp, c("GENE"))
  gene_exp_t <- t(gene_exp[shared_genes, all_samples, with=FALSE])
  setattr(gene_exp_t, 'dimnames', list(NULL, shared_genes))

  exp_grouptest_driven_by_cna <- function(cna_type = 'gain') {
    if (cna_type == 'gain') {
      cna_driven_genes <- gol_ratio_table[Gain > gain_ratio, .(GENE)][[1]]
    } else if (cna_type == 'loss') {
      cna_driven_genes <- gol_ratio_table[Loss > loss_ratio, .(GENE)][[1]]
    } else {
      stop("Unknown cna_type, should be either 'gain' or 'loss'")
    }
    cna_driven_exp <- gene_exp_t[, c(cna_driven_genes)]
    cna_driven_cna <- gene_cna_t[, c(cna_driven_genes)]

    gain_mask <- cna_driven_cna == 1
    normal_mask <- cna_driven_cna == 0
    loss_mask <- cna_driven_cna == -1

    num_samples <- length(all_samples)

    dt <- as.data.table(adply(
      cna_driven_genes,
      1,
      function(gene) {
        cat("Gene:", gene)
        g_exp <- cna_driven_exp[, c(gene)]
        g_gain_mask <- gain_mask[, c(gene)]
        g_normal_mask <- normal_mask[, c(gene)]
        g_loss_mask <- loss_mask[, c(gene)]

        g_gain_exp <- g_exp[g_gain_mask]
        g_normal_exp <- g_exp[g_normal_mask]
        g_loss_exp <- g_exp[g_loss_mask]

        if (cna_type == 'gain') {
          if (length(g_gain_exp) == num_samples) {
            p_value <- NA
          } else {
            p_value = t.test(g_gain_exp, c(g_normal_exp, g_loss_exp))$p.value
          }
        } else {
          if (length(g_loss_exp) == num_samples) {
            p_value <- NA
          } else {
            p_value = t.test(g_loss_exp, c(g_normal_exp, g_gain_exp))$p.value
          }
        }

        ret_val <- data.table(
          GENE = gene,
          p_value = p_value,
          gol_ratio_table[gene, !"GENE", with=FALSE],
          gain_exp_mean = mean(g_gain_exp, na.rm = TRUE),
          normal_exp_mean = mean(g_normal_exp, na.rm = TRUE),
          loss_exp_mean = mean(g_loss_exp, na.rm = TRUE)
        )
      },
      .id = NULL,
      .progress = 'time'
    ))
    dt[, fdr := p.adjust(p_value, method = "fdr")]
  }

  dt_gain <- exp_grouptest_driven_by_cna(cna_type = "gain")
  dt_loss <- exp_grouptest_driven_by_cna(cna_type = "loss")

  return (list(
    dt_gain = dt_gain,
    dt_loss = dt_loss
  ))
}
