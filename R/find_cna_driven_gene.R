#' Perform an integrated analysis of gene expression (GE) and copy number
#' alteration (CNA)
#'
#' The function finds CNA-driven differentially expressed gene and returns the
#' corresponding p-value, false discovery rate, and associated statistics. The
#' result includes three tables which collects information for gain-, loss-, and
#' both-driven genes.
#'
#' The gene is considered CNA-gain if the proportion of the sample exhibiting
#' gain exceeds the threshold \code{gain_prop}, that is, number of samples
#' having \code{gain_loss} = 1. Reversely, the gene is considered CNA-loss if
#' \%samples that \code{gain_loss} = -1 is below a given threshold
#' \code{loss_prop}.
#'
#' When performing the t-test, sample grouping depends on the analysis scenario
#' being either CNA-gain or CNA-loss driven. In CNA-gain driven scenario, two
#' groups, CNA-gain and the other samples, are made. In CNA-loss driven
#' scenario, group CNA-loss and the others are made. Genes that appear in both
#' scenarios will be collected into a third table and excluded from their
#' original tables.
#'
#' See the vignette for usage of this function by a thorough example.
#'
#' @param gene_cna Joint CNA table from \link{create_gene_cna}.
#' @param gene_exp Joint gene expression table from \link{create_gene_exp}.
#' @param gain_prop Minimum proportion of the gain samples to be consider
#'   CNA-gain. Default is 0.2.
#' @param loss_prop Minimum proportion of the loss samples to be consider
#'   CNA-loss. Default is 0.2.
#'
#' @param progress Whether to display a progress bar. By default \code{TRUE}.
#' @param progress_width The text width of the shown progress bar. By default is
#'   48 chars wide.
#' @param parallel Enable parallelism by plyr. One has to specify a parallel
#'   engine beforehand. See example for more information.
#'
#' @return List of three data.table objects for CNA-driven scenarios: gain,
#'   loss, and both, which can be accessed by names: `gain_driven`,
#'   `loss_driven` and `both`.
#'
#' @examples
#' require(data.table)
#'
#' ## Create gene_exp and gene_cna manually. The following shows an example
#' ## consisting of 3 genes (BRCA2, TP53, and GNPAT) and 5 samples (A to E).
#'
#' gene_exp <- data.table(
#'     GENE = c("BRCA2", "TP53", "GNPAT"),
#'     A = c(-0.95, 0.89, 0.21), B = c(1.72, -0.05, NA),
#'     C = c(-1.18, 1.15, 2.47), D = c(-1.24, -0.07, 1.2),
#'     E = c(1.01, 0.93, 1.54)
#' )
#' gene_cna <- data.table(
#'     GENE = c("BRCA2", "TP53", "GNPAT"),
#'     A = c(1, 1, NA), B = c(-1, -1, 1),
#'     C = c(1, -1, 1), D = c(1, -1, -1),
#'     E = c(0, 0, -1)
#' )
#'
#'
#' ## Find CNA-driven genes
#'
#' cna_driven_genes <- find_cna_driven_gene(
#'     gene_cna, gene_exp, progress=FALSE
#' )
#'
#' # Gain driven genes
#' cna_driven_genes$gain_driven
#'
#' # Loss driven genes
#' cna_driven_genes$loss_driven
#'
#' # Gene shown in both gain and loss records
#' cna_driven_genes$both
#'
#'
#' @import data.table
#' @import plyr
#' @export
find_cna_driven_gene <- function(
    gene_cna, gene_exp,
    gain_prop = 0.2, loss_prop = 0.2,
    progress = TRUE, progress_width = 32,
    parallel = FALSE
) {
    all_samples <- colnames(gene_cna)[-1]
    # get shared genes
    shared_genes <- intersect(gene_cna[, .(GENE)][[1]], gene_exp[, .(GENE)][[1]])

    # transpose so columns to be gene-wise, easier to compute
    setkeyv(gene_cna, c("GENE"))
    gene_cna_t <- t(gene_cna[shared_genes, all_samples, with=FALSE])
    # see data.table's setattr doc page.
    # The following equals to
    #     dinames(gene_cna_t) <- list(NULL, shared_genes)
    # setattr prevent from creating unneccessary object copies.
    setattr(gene_cna_t, 'dimnames', list(NULL, shared_genes))

    if (progress) message("Computing gain/loss sample proportion ...\n")
    num_sample <- length(all_samples)
    gol_prop_table <- data.table(
        Gain = colSums(gene_cna_t == 1, na.rm = TRUE) / num_sample,
        Loss = colSums(gene_cna_t == -1, na.rm = TRUE) / num_sample,
        Normal = colSums(gene_cna_t == 0, na.rm = TRUE) / num_sample
    )
    gol_prop_table[, GENE:=shared_genes]
    setkeyv(gol_prop_table, c("GENE"))

    setkeyv(gene_exp, c("GENE"))
    gene_exp_t <- t(gene_exp[shared_genes, all_samples, with=FALSE])
    setattr(gene_exp_t, 'dimnames', list(NULL, shared_genes))

    # define a subroutine to re-use same part computing gain/loss driven genes
    exp_grouptest_driven_by_cna <- function(cna_type = 'gain') {
        if (progress) {
            progress_bar <- progress_text(width = progress_width)
        } else {
            progress_bar <- "none"
        }

        if (cna_type == 'gain') {
            cna_driven_genes <- gol_prop_table[Gain > gain_prop, .(GENE)][[1]]
        } else if (cna_type == 'loss') {
            cna_driven_genes <- gol_prop_table[Loss > loss_prop, .(GENE)][[1]]
        } else {
            stop("Unknown cna_type, should be either 'gain' or 'loss'")
        }
        cna_driven_exp <- gene_exp_t[, c(cna_driven_genes), drop=FALSE]
        cna_driven_cna <- gene_cna_t[, c(cna_driven_genes), drop=FALSE]

        gain_mask <- cna_driven_cna == 1
        normal_mask <- cna_driven_cna == 0
        loss_mask <- cna_driven_cna == -1

        num_samples <- length(all_samples)

        # compute the p-value between gain(loss) vs rest and their avg. gene expression
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
                    if (length(g_gain_exp) <= 1 ||
                        length(c(g_normal_exp, g_loss_exp)) <= 1) {
                        p_value <- NA
                    } else {
                        p_value <- t.test(g_gain_exp, c(g_normal_exp, g_loss_exp))$p.value
                    }
                    vs_rest_exp_diff <- mean(g_gain_exp) - mean(c(g_normal_exp, g_loss_exp))
                } else {
                    if (length(g_loss_exp) <= 1 ||
                        length(c(g_normal_exp, g_gain_exp)) <= 1) {
                        p_value <- NA
                    } else {
                        p_value <- t.test(g_loss_exp, c(g_normal_exp, g_gain_exp))$p.value
                    }
                    vs_rest_exp_diff <- mean(g_loss_exp) - mean(c(g_normal_exp, g_gain_exp))
                }

                ret_val <- data.table(
                    GENE = gene,
                    p_value = p_value,
                    gol_prop_table[gene, !"GENE", with=FALSE],
                    gain_exp_mean = mean(g_gain_exp),
                    normal_exp_mean = mean(g_normal_exp),
                    loss_exp_mean = mean(g_loss_exp),
                    vs_rest_exp_diff = vs_rest_exp_diff
                )
            },
            .id = NULL,
            .progress = progress_bar,
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
            c("gain_sample_prop", "normal_sample_prop", "loss_sample_prop")
        )
        return(cna_driven_dt)
    }

    # use the subroutine
    if (progress) message("Computing CNA gain driven gene records ... \n")
    dt_gain <- exp_grouptest_driven_by_cna(cna_type = "gain")

    if (progress) message("Computing CNA loss driven gene records ... \n")
    dt_loss <- exp_grouptest_driven_by_cna(cna_type = "loss")

    # [.data.table requires key columns
    setkey(dt_gain, "GENE")
    setkey(dt_loss, "GENE")

    # take out genes shown in both gain and loss table
    dt_both <- dt_gain[, .(GENE, p_value, fdr, vs_rest_exp_diff)][dt_loss, nomatch=0]
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
