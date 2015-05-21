CNAtoGene <- function(file.pattern, directory, tcga = TRUE, cna.gain.threshold, 
    cna.loss.threshold, col.sample, col.chrmosome, col.startloci, col.endloci, 
    col.seg.mean, outputList = FALSE) {
    if (missing(directory)) {
        path_cna <- getwd()
    } else {
        path_cna <- directory
    }
    
    load("C:/Users/Daphne/Documents/Project/hg19DBNM.rda")
    GeneList <- unique(hg19DBNM[, 6])
    iniV = c(rep(0, length(GeneList)))
    cna2geneList <- vector()
    CNAwholeGeneList <- vector()
    
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
            files_cna <- list.files(path = path_cna, pattern = paste("*\\", file.pattern, 
                sep = ""))
        }
    }
    
    if (missing(col.sample)) {
        checkCOLloc <- scan(paste(path_cna, files_cna[1], sep = "/"), "", nlines = 1)
        idx_colsample <- grep("sample", checkCOLloc, ignore.case = T)
        idx_colchr <- grep("chr", checkCOLloc, ignore.case = T)
        idx_colstartloci <- grep("start", checkCOLloc, ignore.case = T)
        idx_colendloci <- grep("end", checkCOLloc, ignore.case = T)
        idx_colsegmean <- grep("mean", checkCOLloc, ignore.case = T)
    } else {
        idx_colsample <- col.sample
        idx_colchr <- col.chrmosome
        idx_colstartloci <- col.startloci
        idx_colendloci <- col.endloci
        idx_colsegmean <- col.seg.mean
    }
    
    if (missing(cna.gain.threshold)) {
        cna_gain_thres <- log(2.5, 2) - 1
    } else {
        cna_gain_thres <- cna.gain.threshold
    }
    if (missing(cna.loss.threshold)) {
        cna_loss_thres <- log(1.5, 2) - 1
    } else {
        cna_loss_thres <- cna.loss.threshold
    }
    
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
    CNAtoGeneList <- CNAwholeGeneList
    save(CNAtoGeneList, file = paste(path_cna, "CNAtoGeneList.rda", sep = "/"))
    
    if (outputList) {
        write.csv(CNAtoGeneList, file = paste(path_cna, "CNAtoGeneList.csv", sep = "/"), 
            row.names = F)
    }
} 
