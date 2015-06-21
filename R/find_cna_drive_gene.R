#' @import data.table
#' @import plyr
#' @export
find_cna_driven_gene <- function(
  gene_cna, gene_exp,
  gain_ratio = 0.2, loss_ratio = 0.2,
  progress = TRUE, progress_width = 32,
  parallel = FALSE
) {
  all_samples <- colnames(gene_cna)[-1]
  # get shared genes
  shared_genes <- intersect(
    gene_cna[, .(GENE)][[1]],
    gene_exp[, .(GENE)][[1]]
  )

  # transpose so columns to be gene-wise, easier to compute
  setkeyv(gene_cna, c("GENE"))
  gene_cna_t <- t(gene_cna[shared_genes, all_samples, with=FALSE])
  # see data.table's setattr doc page.
  # The following equals to
  #     dinames(gene_cna_t) <- list(NULL, shared_genes)
  # setattr prevent from creating unneccessary object copies.
  setattr(gene_cna_t, 'dimnames', list(NULL, shared_genes))

  cat("Computing gain/loss sample ratio ...\n")
  num_sample <- length(all_samples)
  gol_ratio_table <- data.table(
    Gain = colSums(gene_cna_t == 1) / num_sample,
    Loss = colSums(gene_cna_t == -1) / num_sample,
    Normal = colSums(gene_cna_t == 0) / num_sample
  )
  # gol_ratio_table <- as.data.table(aaply(
  #   gene_cna_t,
  #   2,
  #   function(one_gene) {
  #     c(Gain = sum(one_gene == 1),
  #       Loss = sum(one_gene == -1),
  #       Normal = sum(one_gene == 0))},
  #   .progress = progress_text(width=progress_width)
  # ) / length(all_samples))
  gol_ratio_table[, GENE:=shared_genes]
  setkeyv(gol_ratio_table, c("GENE"))

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
    cna_driven_exp <- gene_exp_t[, c(cna_driven_genes), drop=FALSE]
    cna_driven_cna <- gene_cna_t[, c(cna_driven_genes), drop=FALSE]

    gain_mask <- cna_driven_cna == 1
    normal_mask <- cna_driven_cna == 0
    loss_mask <- cna_driven_cna == -1

    num_samples <- length(all_samples)

    cna_driven_dt <- as.data.table(adply(
      cna_driven_genes,
      1,
      function(gene) {
        g_exp <- cna_driven_exp[, c(gene)]
        g_gain_mask <- gain_mask[, c(gene)]
        g_normal_mask <- normal_mask[, c(gene)]
        g_loss_mask <- loss_mask[, c(gene)]

        g_gain_exp <- na.omit(g_exp[g_gain_mask])
        g_normal_exp <- na.omit(g_exp[g_normal_mask])
        g_loss_exp <- na.omit(g_exp[g_loss_mask])

        if (cna_type == 'gain') {
          if (length(g_gain_exp) <= 1 || length(c(g_normal_exp, g_loss_exp)) <= 1) {
            p_value <- NA
            vs_rest_exp_diff <- NA
          } else {
            p_value = t.test(g_gain_exp, c(g_normal_exp, g_loss_exp))$p.value
            vs_rest_exp_diff <- mean(g_gain_exp) - mean(c(g_normal_exp, g_loss_exp))
          }
        } else {
          if (length(g_loss_exp) <= 1 || length(c(g_normal_exp, g_gain_exp)) <= 1) {
            p_value <- NA
            vs_rest_exp_diff <- NA
          } else {
            p_value = t.test(g_loss_exp, c(g_normal_exp, g_gain_exp))$p.value
            vs_rest_exp_diff <- mean(g_loss_exp) - mean(c(g_normal_exp, g_gain_exp))
          }
        }

        ret_val <- data.table(
          GENE = gene,
          p_value = p_value,
          gol_ratio_table[gene, !"GENE", with=FALSE],
          gain_exp_mean = mean(g_gain_exp),
          normal_exp_mean = mean(g_normal_exp),
          loss_exp_mean = mean(g_loss_exp),
          vs_rest_exp_diff = vs_rest_exp_diff
        )
      },
      .id = NULL,
      .progress = progress_text(width=progress_width),
      .parallel = parallel
    ))
    fdr_adjusted_p <- p.adjust(cna_driven_dt$p_value, method = "fdr")
    set(cna_driven_dt, NULL, "fdr", fdr_adjusted_p)
    setcolorder(
      cna_driven_dt,
      c("GENE", "p_value", "fdr",
        "Gain", "Normal", "Loss",
        "gain_exp_mean", "normal_exp_mean", "loss_exp_mean",
        "vs_rest_exp_diff")
    )
    setnames(
      cna_driven_dt,
      c("Gain", "Normal", "Loss"),
      c("gain_sample_ratio", "normal_sample_ratio", "loss_sample_ratio")
    )
    return(cna_driven_dt)
  }

  cat("Computing CNA gain driven gene records ... \n")
  dt_gain <- exp_grouptest_driven_by_cna(cna_type = "gain")
  cat("Computing CNA loss driven gene records ... \n")
  dt_loss <- exp_grouptest_driven_by_cna(cna_type = "loss")

  # [.data.table requires key columns
  setkey(dt_gain, "GENE")
  setkey(dt_loss, "GENE")

  dt_both <- dt_gain[, .(GENE, p_value, fdr, vs_rest_exp_diff)][dt_loss, nomatch=0]
  # dt_both <- merge(dt_gain[, .(GENE, p_value, fdr)], dt_loss, by="GENE", all=TRUE)
  setnames(
    dt_both,
    c("p_value", "fdr", "vs_rest_exp_diff",
      "i.p_value", "i.fdr", "i.vs_rest_exp_diff"),
    c("gain_p_value", "gain_fdr", "gain_vs_rest_exp_diff",
      "loss_p_value", "loss_fdr", "loss_vs_exp_diff")
  )
  # move *_vs_rest_exp_diff to the end
  setcolorder(
    dt_both,
    c(colnames(dt_both)[-c(4, 13)], colnames(dt_both)[c(4, 13)])
  )
  both_driven_genes <- dt_both[, .(GENE)]

  dt_gain <- dt_gain[!both_driven_genes]
  dt_loss <- dt_loss[!both_driven_genes]

  # unset keys
  setkey(dt_gain, NULL)
  setkey(dt_loss, NULL)
  setkey(dt_both, NULL)
  # sorted by ascending fdr
  setorderv(dt_gain, "fdr", 1, na.last=TRUE)
  setorderv(dt_loss, "fdr", 1, na.last=TRUE)
  setorderv(dt_both, c("gain_fdr", "loss_fdr"), c(1, 1), na.last=TRUE)

  return (list(
    gain_driven = dt_gain,
    loss_driven = dt_loss,
    both = dt_both
  ))
}
