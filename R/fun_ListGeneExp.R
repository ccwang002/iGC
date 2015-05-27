#' Create an overall gene expression list of entire input samples
#'
#' \code{ListGeneExp} returns a gene expression list of all input files from The Cancer Genome Atlas (TCGA),
#' Gene Expression Omnibus (GEO), and other sources.
#'
#' @section geneExpList function:
#' The ListGeneExp function is used to rearrange the input gene expression files to a gene expression
#' list of entire samples.
#'
#' @param file
#'        Filename for the gene list of input sample.
#' @param directory
#'        The directory gene list of input sample locates.
#' @param tcga
#'        True if TCGA's data is used. \strong{Default value is TRUE}.
#' @param probe
#'        True if \strong{the gene list of input sample has probe information}. \strong{Default value is FALSE}.
#' @param outputList
#'        True if \emph{user wants to get the rearranged output list of all input files in csv format}.
#'        \strong{Default value is TRUE}.
#' @return  gene list.
#' @examples
#' \dontrun{
#' ListGeneExp('my_gene_list.txt', "C:/Users/Documents/Project/Input_samples", True, False, True)
#' }
#' @seealso
#' iGC-an integrated analysis package of Gene expression and Copy number alteration
#' @export
ListGeneExp <- function(
  file, directory,
  tcga = TRUE,
  probe = FALSE,
  outputList = TRUE
) {
    if (missing(directory)) {
        path <- getwd()
    } else {
        path <- directory
    }

    GeneExp <- vector()
    if (!tcga) {
        GeneExp <- read.csv(paste(path, file, sep = "/"), header = TRUE)
        if (probe) {
            colnames(GeneExp[1]) <- c("Probe")
            colnames(GeneExp[2]) <- c("Gene")
        } else if (!probe) {
            colnames(GeneExp[1]) <- c("Gene")
        }
    } else if (tcga || missing(file)) {
        files <- list.files(path = path, pattern = "*.data")
        for (file in files) {
            bindtemp <- vector()
            title <- vector()
            title <- scan(paste(path, file, sep = "/"), "", nlines = 1)
            bindtemp <- read.table(paste(path, file, sep = "/"), skip = 2)
            if (sapply(bindtemp, class)[2] == "factor") {
                bindtemp[, 2] <- as.numeric(as.character(bindtemp[, 2]))
            }
            names(bindtemp) <- c("GENE", title[3])
            if (length(GeneExp) == 0) {
                GeneExp <- bindtemp
            } else {
                GeneExp <- merge(GeneExp, bindtemp, by = "GENE", sort = F)
            }
        }
    }

    save(GeneExp, file = paste(path, "GeneExp.rda", sep = "/"))

    if (outputList) {
        write.csv(GeneExp, file = paste(path, "GeneExp.csv", sep = "/"), row.names = F)
    }
}
