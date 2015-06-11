#' Map CNA data to human genes using human genome version 19 as reference
#'
#' \code{CNAtoGene} returns a gene list of its CNA status, CNA-gain (G), CNA-loss (L) and neutral (N),
#' in all input samples and the ratio of the number of CN changed samples (G or L) to the number of total samples.
#'
#'
#' @section CNAtoGene function:
#' The CNAtoGene function maps CNA data to human genes and then defines the mapped human genes as CN gain or loss
#' based on the CN threshold, whose default values are set as 2.5 for gain and 1.5 for loss.
#' These mapped genes will be assigned values in +1, -1 or 0, where +1 stands for CNA-gain, -1 stands for CNA-loss
#' and 0 stands for neutral.
#'
#' @param file.pattern
#'        Filename or file pattern of the input CNA samples.
#' @param directory
#'        The directory list of input CNA sample locates.
#' @param tcga
#'        True if TCGA's data is used. \strong{Default value is TRUE}.
#' @param cna.gain.threshold
#'        CNA expression threshold for defining as CNA-gain, where \strong{default value is set as 2.5 for gain}.
#' @param cna.loss.threshold
#'        CNA expression threshold for defining as CNA-loss, where \strong{default value is set as 1.5 for loss}.
#' @param col.sample
#'        The column number of sample.
#' @param col.chromosome
#'        The column number of chromosome.
#' @param col.startloci
#'        The column number of starting position of each copy number segmentation.
#' @param col.endloci
#'        The column number of ending position of each copy number segmentation.
#' @param col.seg.mean
#'        The column number of CNA mean of each copy number segmentation.
#' @param outputList
#'        True if \emph{user wants to get the output list of mapped genes with CNA status in csv format}.
#'        \strong{Default value is TRUE}.
#' @return  mapped gene list with CNA status
#' @examples
#' \dontrun{
#' CNAtoGene('*\\.hg19.seg.txt', "C:/Users/Documents/Project/Input_samples", True, 2.7, 1.3, 1, 2, 3, 4, 6, True)
#' }
#' @seealso
#' iGC-an integrated analysis package of Gene expression and Copy number alteration
#' @export
CNAtoGene <- function(file.pattern, directory, tcga = TRUE, cna.gain.threshold,
    cna.loss.threshold, col.sample, col.chromosome, col.startloci, col.endloci,
    col.seg.mean, outputList = True) {
    if (missing(directory)) {
        path_cna <- getwd()
    } else {
        path_cna <- directory
    }
    print(path_cna)

    data("hg19DBNM", package = "iGC", envir = environment())
    GeneList <- unique(hg19DBNM[, 6])
    iniV = c(rep(0, length(GeneList)))

    if (!tcga) {
        if (missing(file.pattern)) {
            files_cna <- list.files(path = path_cna)
        } else {
            files_cna <- list.files(path = path_cna, file.pattern)
        }
    } else if (tcga) {
        if (missing(file.pattern)) {
            files_cna <- list.files(path = path_cna, pattern = "*\\.hg19.seg.txt")
            if (length(files_cna) == 0) {
                files_cna <- list.files(path = path_cna)
            }
        } else {
            files_cna <- list.files(path = path_cna, pattern = file.pattern)
        }
    }
    print(files_cna)

    if (missing(col.sample)) {
        checkCOLloc <- scan(paste(path_cna, files_cna[1], sep = "/"), "", nlines = 1)
        idx_colsample <- grep("sample", checkCOLloc, ignore.case = T)
        idx_colchr <- grep("chr", checkCOLloc, ignore.case = T)
        idx_colstartloci <- grep("start", checkCOLloc, ignore.case = T)
        idx_colendloci <- grep("end", checkCOLloc, ignore.case = T)
        idx_colsegmean <- grep("mean", checkCOLloc, ignore.case = T)
    } else {
        idx_colsample <- col.sample
        idx_colchr <- col.chromosome
        idx_colstartloci <- col.startloci
        idx_colendloci <- col.endloci
        idx_colsegmean <- col.seg.mean
    }

    if (missing(cna.gain.threshold)) {
        cna_gain_thres <- log(2.5, 2) - 1
    } else {
        cna_gain_thres <- cna.gain.threshold
    }
    print(cna_gain_thres)
    if (missing(cna.loss.threshold)) {
        cna_loss_thres <- log(1.5, 2) - 1
    } else {
        cna_loss_thres <- cna.loss.threshold
    }
    print(cna_loss_thres)

    CNAwholeGeneList <- vector()
    for (file in files_cna) {
        cna_main <- vector()
        title_cna <- vector()
        cna2gene <- vector()
        wholeGeneList <- data.frame(GeneList, iniV)
        if (tcga) {
            title_cna <- scan(paste(path_cna, file, sep = "/"), "", skip = 1, nlines = 1)
            cna_main <- read.table(paste(path_cna, file, sep = "/"), header = T)
        } else {
            idx.tit <- regexpr(".csv", file)
            title_cna <- substr(file, 1, idx.tit[1] - 1)
            cna_main <- read.csv(paste(path_cna, file, sep = "/"), header = T)
        }
        for (ii in 1:dim(cna_main)[1]) {
            idx_cn <- vector()
            idx_db <- vector()
            idx_cn = cna_main[ii, idx_colchr]
            if (idx_cn == 23 || idx_cn == "X") {
                idx_db = which(hg19DBNM[, 2] == "X")
            } else if (idx_cn == 24 || idx_cn == "Y") {
                idx_db = which(hg19DBNM[, 2] == "Y")
            } else idx_db = which(hg19DBNM[, 2] == as.character(idx_cn))


            if (length(idx_db) > 0) {
                idx_srt = head(which(cna_main[ii, idx_colstartloci] < hg19DBNM[idx_db[1]:tail(idx_db,
                  1), 4]), 1)
                idx_end = tail(which(cna_main[ii, idx_colendloci] > hg19DBNM[idx_db[1]:tail(idx_db,
                  1), 3]), 1)

                if (!is.na(idx_srt > 0 && idx_end > 0) && idx_end > idx_srt) {
                  GL = unique(hg19DBNM[idx_db[idx_srt:idx_end], 6])

                  if (is.na(cna_main[ii, idx_colsegmean])) {
                    gol = 0
                  } else if (cna_main[ii, idx_colsegmean] > cna_gain_thres) {
                    gol = 1
                  } else if (cna_main[ii, idx_colsegmean] < cna_loss_thres) {
                    gol = -1
                  } else gol = 0

                  gol_all = c(rep(gol, length(GL)))

                  cna2gene_temp <- data.frame(GL, gol_all)
                  names(cna2gene_temp) <- c("Genes", title_cna[1])

                  cna2gene <- rbind(cna2gene, cna2gene_temp)
                }
            }
        }
        idx_GL = which(cna2gene[, 2] != 0)
        idx_map = match(cna2gene[idx_GL, 1], wholeGeneList[, 1])
        wholeGeneList[idx_map, 2] <- cna2gene[idx_GL, 2]
        names(wholeGeneList) <- c("GeneList", title_cna[1])
        if (length(CNAwholeGeneList) == 0) {
            CNAwholeGeneList <- wholeGeneList
        } else {
            CNAwholeGeneList <- merge(CNAwholeGeneList, wholeGeneList, by = "GeneList",
                sort = F)
        }

    }
    CNAtoGeneList <- CNAwholeGeneList

    path.out <- getwd()
    print(path.out)
    save(CNAtoGeneList, file = paste(path.out, "CNAtoGeneList.rda", sep = "/"))

    if (outputList) {
        write.csv(CNAtoGeneList, file = paste(path.out, "CNAtoGeneList.csv", sep = "/"),
            row.names = F)
    }
}
