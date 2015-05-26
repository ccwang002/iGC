#' The human genome reference used here is RefSeq transcripts in version hg19, hg19-RefSeq,
#' from UCSC Genome Browser.
#'
#' This reference provides region information, including chromosome number, starting position,
#' ending position, strand and gene symbols, for converting copy number alteration data into human genes.
#'
#' The transcripts with NM marker ID, which are protein-codeing, were selected to be our reference database
#' and provided as hg19DBNM.rda.
#'
#' @format A data frame with 39997 rows and 7 variables:
#' \describe{
#'   \item{marker_id}{RefSeq name with its corrsponding gene symbol}
#'   \item{chromosome}{1-21, X and Y}
#'   \item{start}{starting position, in basepair number}
#'   \item{stop}{ending position, in basepair number}
#'   \item{strand}{positive or negative strand, in + or - symbols}
#'   \item{Gene.Symbol}{Gene name}
#'   \item{Transcript}{RefSeq name}
#' }
#' @source \url{UCSC Genome Browser, http://hgdownload.cse.ucsc.edu/downloads.html}
"hg19DBNM"
