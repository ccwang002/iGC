#' hg19-RefSeq
#'
#' The human genome reference used here is RefSeq transcripts in version hg19
#' from UCSC Genome Browser. The transcripts with NM marker ID, which are
#' protein-codeing, were selected to be our reference database and provided as
#' hg19DBNM.rda.
#'
#' This reference provides region information, including chromosome number,
#' starting position, ending position, strand and gene symbols, for converting
#' copy number alteration data into human genes.
#'
#' @format A data frame with 39997 rows and 7 variables:
#' \describe{
#'   \item{Marker.ID}{RefSeq name with its corrsponding gene symbol}
#'   \item{Chromosome}{1-22, X and Y}
#'   \item{Start}{starting position, in basepair number}
#'   \item{Stop}{ending position, in basepair number}
#'   \item{Strand}{positive or negative strand, in + or - symbols}
#'   \item{Gene.Symbol}{Gene name}
#'   \item{Transcript}{RefSeq name}
#' }
#' @return data.table
#' @source UCSC Genome Browser:
#'   \url{http://hgdownload.cse.ucsc.edu/downloads.html}
"hg19DBNM"
