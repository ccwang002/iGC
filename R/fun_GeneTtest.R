CNAdrivenGene <- function(CNAtoGeneList, GeneExp, probe = FALSE, sample.mapping.file, directory, map.loc.start, 
    map.loc.end, gain.ratio, loss.ratio, tcga = TRUE, outputList = TRUE, output_GSEA = FALSE) {
    if (missing(directory)) {
        path_map <- getwd()
    } else {
        path_map <- directory
    }
    
    CNAwholeGeneList <- read.csv(CNAtoGeneList, header = T, check.names = FALSE)
    
    numGL <- CNAwholeGeneList[, 1]
    iniV = c(rep(0, length(numGL)))
    ListGL <- data.frame(numGL, iniV, iniV, iniV)
    names(ListGL) <- c("GeneList", "numGain", "numLoss", "numNormal")
    
    if (!tcga) {
        if (missing(sample.mapping.file)) {
            CNAwholeGeneList_tumor <- CNAwholeGeneList
        } else {
            files_map <- list.files(path = path_map, pattern = sample.mapping.file)
            cna_map <- vector()
            cna_map <- read.csv(paste(path_map, files_map, sep = "/"), header = T)
            sampleName1 <- colnames(CNAwholeGeneList)
            idx_match <- match(sampleName1, cna_map[, 1])
            sampleName2 <- cna_map[idx_match, 2]
            sampleName <- rbind(sampleName1, as.character(sampleName2))
            sampleName.nona <- sampleName[, complete.cases(sampleName[2, ])]
            idx.nona <- match(sampleName.nona[1, ], sampleName1)
            CNAwholeGeneList_nona <- CNAwholeGeneList[, c(1, idx.nona)]
            colnames(CNAwholeGeneList_nona) <- c("GeneList", paste(as.character(sampleName.nona[2, ])))
            CNAwholeGeneList_tumor <- CNAwholeGeneList_nona
        }
    } else if (tcga) {
        if (missing(sample.mapping.file)) {
            files_map <- list.files(path = path_map, pattern = "FILE_SAMPLE_MAP.txt")
        } else {
            files_map <- list.files(path = path_map, pattern = sample.mapping.file)
        }
        cna_map_input <- vector()
        cna_map_input <- read.table(paste(path_map, files_map, sep = "/"), header = T)
        idx_map <- grep("*\\.hg19.seg.txt", cna_map_input[, 1])
        if (length(idx_map) > 0) {
            cna_map <- cna_map_input[idx_map, ]
        } else {
            cna_map <- cna_map_input
        }
        sampleName1 <- colnames(CNAwholeGeneList)
        idx_match <- match(substr(sampleName1, start = map.loc.start, stop = map.loc.end), substr(cna_map[, 
            1], start = map.loc.start, stop = map.loc.end))
        sampleName2 <- cna_map[idx_match, 2]
        sampleName <- rbind(sampleName1, as.character(sampleName2))
        CNAwholeGeneList_name <- CNAwholeGeneList
        colnames(CNAwholeGeneList_name) <- c("GeneList", paste(as.character(sampleName2[2:length(sampleName2)])))
        idx_tumor <- which(substr(colnames(CNAwholeGeneList_name), start = 14, stop = 14) == "0")
        CNAwholeGeneList_tumor <- CNAwholeGeneList_name[, c(1, idx_tumor)]
    }
    
    for (kk in 1:nrow(CNAwholeGeneList_tumor)) {
        idx_gain <- which(CNAwholeGeneList_tumor[kk, ] == 1)
        idx_loss <- which(CNAwholeGeneList_tumor[kk, ] == -1)
        idx_normal <- which(CNAwholeGeneList_tumor[kk, ] == 0)
        ListGL[kk, 2] <- length(idx_gain)
        ListGL[kk, 3] <- length(idx_loss)
        ListGL[kk, 4] <- length(idx_normal)
    }
    ListGLratio <- data.frame(numGL, ListGL[, 2:4]/(ncol(CNAwholeGeneList_tumor) - 1))
    names(ListGLratio) <- c("GeneList", "ratioGain", "ratioLoss", "ratioNormal")
    
    CNAwholeGeneListNum <- c()
    CNAwholeGeneListNum <- merge(CNAwholeGeneList_tumor, ListGLratio, by = "GeneList", sort = F)
    
    if (missing(gain.ratio)) {
        gain_ratio <- 0.2
    } else {
        gain_ratio <- gain.ratio
    }
    if (missing(loss.ratio)) {
        loss_ratio <- 0.2
    } else {
        loss_ratio <- loss.ratio
    }
    
    idx_thG <- which(CNAwholeGeneListNum[, ncol(CNAwholeGeneListNum) - 2] > gain_ratio)
    idx_thL <- which(CNAwholeGeneListNum[, ncol(CNAwholeGeneListNum) - 1] > loss_ratio)
    
    
    GeneExp <- read.csv(GeneExp, header = T, check.names = FALSE)
    
    if (probe) {
        GeneExp <- GeneExp[, 2:dim(GeneExp)[2]]
    }
    
    thGlist <- CNAwholeGeneListNum[idx_thG, ]
    thLlist <- CNAwholeGeneListNum[idx_thL, ]
    
    idx_Gmat <- match(thGlist[, 1], GeneExp[, 1])
    thGlist_expmat <- thGlist[which(idx_Gmat != "NA"), ]
    thGlist_expNA <- thGlist[which(is.na(idx_Gmat)), ]
    idx_Lmat <- match(thLlist[, 1], GeneExp[, 1])
    thLlist_expmat <- thLlist[which(idx_Lmat != "NA"), ]
    thLlist_expNA <- thLlist[which(is.na(idx_Lmat)), ]
    
    thGlist_gexp <- thGlist_expmat
    thLlist_gexp <- thLlist_expmat
    
    idx_Ggene <- match(thGlist_gexp[, 1], GeneExp[, 1])
    for (pp in 2:(ncol(thGlist_gexp) - 3)) {
        if (tcga) {
            idx_Gsam <- pmatch(substr(colnames(thGlist_gexp[pp]), start = 1, stop = 16), colnames(GeneExp))
        } else {
            idx_Gsam <- match(colnames(thGlist_gexp[pp]), colnames(GeneExp))
        }
        thGlist_gexp[, pp] <- GeneExp[idx_Ggene, idx_Gsam]
    }
    
    idx_Lgene <- match(thLlist_gexp[, 1], GeneExp[, 1])
    for (qq in 2:(ncol(thLlist_gexp) - 3)) {
        if (tcga) {
            idx_Lsam <- pmatch(substr(colnames(thLlist_gexp[qq]), start = 1, stop = 16), colnames(GeneExp))
        } else {
            idx_Lsam <- match(colnames(thLlist_gexp[qq]), colnames(GeneExp))
        }
        thLlist_gexp[, qq] <- GeneExp[idx_Lgene, idx_Lsam]
    }
    
    GeneList_Gain <- thGlist_gexp[, 1]
    pValue_Gain <- c(rep(0, nrow(thGlist_gexp)))
    meanGain <- c(rep(0, nrow(thGlist_gexp)))
    meanLoss <- c(rep(0, nrow(thGlist_gexp)))
    meanNormal <- c(rep(0, nrow(thGlist_gexp)))
    meanDiff_GnonG <- c(rep(0, nrow(thGlist_gexp)))
    
    for (rr in 1:nrow(thGlist_gexp)) {
        idx_GG <- which(thGlist_expmat[rr, 1:(ncol(thGlist) - 3)] == 1)
        idx_GN <- which(thGlist_expmat[rr, 1:(ncol(thGlist) - 3)] == 0)
        idx_GL <- which(thGlist_expmat[rr, 1:(ncol(thGlist) - 3)] == -1)
        idx_nonG <- c(idx_GN, idx_GL)
        exp_GG <- thGlist_gexp[rr, idx_GG]
        exp_GN <- thGlist_gexp[rr, idx_GN]
        exp_GL <- thGlist_gexp[rr, idx_GL]
        exp_nonG <- thGlist_gexp[rr, idx_nonG]
        meanGain[rr] <- mean(as.numeric(exp_GG))
        meanNormal[rr] <- mean(as.numeric(exp_GN))
        meanLoss[rr] <- mean(as.numeric(exp_GL))
        meanDiff_GnonG[rr] <- mean(as.numeric(exp_GG)) - mean(as.numeric(exp_nonG))
        if (thGlist_gexp[rr, (ncol(thGlist_gexp) - 2)] == 1) 
            pValue_Gain[rr] <- NA else pValue_Gain[rr] <- t.test(exp_GG, exp_nonG)$p.value
    }
    
    fdr_Gain <- p.adjust(pValue_Gain, method = "fdr")
    Gain_output <- data.frame(GeneList_Gain, meanGain, meanLoss, meanNormal, meanDiff_GnonG, thGlist_expmat[, 
        (ncol(thGlist) - 2):ncol(thGlist)], pValue_Gain, fdr_Gain)
    Gain_sorted <- Gain_output[order(Gain_output$pValue_Gain), ]
    save(Gain_sorted, file = paste(path_map, "Gain_sorted.rda", sep = "/"))
    
    
    GeneList_Loss <- thLlist_gexp[, 1]
    pValue_Loss <- c(rep(0, nrow(thLlist_gexp)))
    meanLoss <- c(rep(0, nrow(thLlist_gexp)))
    meanGain <- c(rep(0, nrow(thLlist_gexp)))
    meanNormal <- c(rep(0, nrow(thLlist_gexp)))
    meanDiff_LnonL <- c(rep(0, nrow(thLlist_gexp)))
    
    
    for (ss in 1:nrow(thLlist_gexp)) {
        idx_LL <- which(thLlist_expmat[ss, 1:(ncol(thLlist) - 3)] == -1)
        idx_LG <- which(thLlist_expmat[ss, 1:(ncol(thLlist) - 3)] == 1)
        idx_LN <- which(thLlist_expmat[ss, 1:(ncol(thLlist) - 3)] == 0)
        idx_nonL <- c(idx_LG, idx_LN)
        exp_LL <- thLlist_gexp[ss, idx_LL]
        exp_LG <- thLlist_gexp[ss, idx_LG]
        exp_LN <- thLlist_gexp[ss, idx_LN]
        exp_nonL <- thLlist_gexp[ss, idx_nonL]
        meanLoss[ss] <- mean(as.numeric(exp_LL))
        meanGain[ss] <- mean(as.numeric(exp_LG))
        meanNormal[ss] <- mean(as.numeric(exp_LN))
        meanDiff_LnonL[ss] <- mean(as.numeric(exp_LL)) - mean(as.numeric(exp_nonL))
        
        if (thLlist_gexp[ss, (ncol(thLlist_gexp) - 1)] == 1) 
            pValue_Loss[ss] <- NA else pValue_Loss[ss] <- t.test(exp_LL, exp_nonL)$p.value
        
    }
    
    fdr_Loss <- p.adjust(pValue_Loss, method = "fdr")
    Loss_output <- data.frame(GeneList_Loss, meanGain, meanLoss, meanNormal, meanDiff_LnonL, thLlist_expmat[, 
        (ncol(thLlist) - 2):ncol(thLlist)], pValue_Loss, fdr_Loss)
    Loss_sorted <- Loss_output[order(Loss_output$pValue_Loss), ]
    save(Loss_sorted, file = paste(path_map, "Loss_sorted.rda", sep = "/"))
    
    
    idx_both <- match(Gain_output[, 1], Loss_output[, 1])
    idx_both_gain <- which(!is.na(idx_both))
    output_both_gain <- Gain_output[idx_both_gain, ]
    idx_both_loss <- match(output_both_gain[, 1], Loss_output[, 1])
    output_both_loss <- Loss_output[idx_both_loss, ]
    Both_GainLoss_sorted <- data.frame(output_both_gain, output_both_loss)
    save(Both_GainLoss_sorted, file = paste(path_map, "Both_GainLoss_sorted.rda", sep = "/"))
    
    
    if (outputList) {
        write.csv(Gain_sorted, file = paste(path_map, "Gain_sorted.csv", sep = "/"), row.names = F)
        write.csv(Loss_sorted, file = paste(path_map, "Loss_sorted.csv", sep = "/"), row.names = F)
        write.csv(Both_GainLoss_sorted, file = paste(path_map, "Both_GainLoss_sorted.csv", sep = "/"), row.names = F)
    }
    
    if (output_GSEA) {
        Gain_sorted_gsea <- Gain_sorted[, 1:5]
        write.table(Gain_sorted_gsea, file = paste(path_map, "Gain_sorted_gsea.txt", sep = "/"), row.names = F, 
            sep = "\t")
        Loss_sorted_gsea <- Loss_sorted[, 1:5]
        write.table(Loss_sorted_gsea, file = paste(path_map, "Loss_sorted_gsea.txt", sep = "/"), row.names = F, 
            sep = "\t")
        Both_GainLoss_sorted_gsea <- data.frame(Both_GainLoss_sorted[, 1:5], Both_GainLoss_sorted[, 12:15])
        write.table(GBoth_GainLoss_sorted_gsea, file = paste(path_map, "Both_GainLoss_sorted_gsea.txt", sep = "/"), 
            row.names = F, sep = "\t")
    }
}












 
