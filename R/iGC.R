#' iGC: an integrated analysis package of gene expression and copy number alteration
#'
#' The iGC package is used to identify CNA-driven differentially expressed genes.
#' The iGC package provides three categories of important functions:
#' ListGeneExp, CNAtoGene and CNAdrivenGene.
#'
#' @section ListGeneExp function:
#' The ListGeneExp function is used to rearrange the input gene expression files to a gene expression
#' list of entire samples.
#'
#' @section CNAtoGene function:
#' The CNAtoGene function maps CNA data to human genes and then defines the mapped human genes as CN gain or loss
#' based on the CN threshold, whose default values are set as 2.5 for gain and 1.5 for loss. These mapped genes will
#' be assigned values in +1, -1 or 0, where +1 stands for CNA-gain, -1 stands for CNA-loss and 0 stands for neutral.
#'
#' @section CNAdrivenGene function:
#' The CNAdrivenGene function identifies CNA-driven differentially expressed genes. The input mapped genes
#' remain for further analyses if its ratio of the number of CN changed samples, CNA-gain (G) or CNA-loss (L),
#' to the number of total samples is larger than a given threshold. Here the default setting is that only genes
#' showing CNAs in at least 20% of the samples will be analyzed further. Then, statistical tests, T-test and
#' Wilcoxon rank sum test, are performed in the GE level by classifying the samples as G and L plus Nertral (N)
#' groups or L and G plus N groups, depending on the CN of the interested gene increases or decreases.
#'
#' @docType package
#' @name iGC
NULL
