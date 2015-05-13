geneExpList <- function(file, directory, tcga=TRUE, probe=FALSE, outputList=FALSE){
	if (missing(directory)){
		path <- getwd()
	}else{
		path <- directory
	}	

	GeneExp <- vector()
	if (!tcga){
		GeneExp <- read.csv(paste(path, file, sep="/"), header=TRUE)
		if (probe){
        		colnames(GeneExp[1]) <- c("Probe")
        		colnames(GeneExp[2]) <- c("Gene")
		}else if (!probe){
        		colnames(GeneExp[1]) <- c("Gene")
      	}
	} else if(tcga || missing(file)){
		files <- list.files(path=path, pattern="*.data")
		for(file in files){
			bindtemp <- vector()
    			title <- vector()
    			title <- scan(paste(path, file, sep="/"), "", nlines=1)
    			bindtemp <- read.table(paste(path, file, sep="/"), skip=2)
    			if (sapply(bindtemp,class)[2]=="factor"){
        			bindtemp[,2] <- as.numeric(as.character(bindtemp[,2]))
			}
			names(bindtemp) <- c("GENE", title[3])
    			if (length(GeneExp)==0){
        			GeneExp <- bindtemp
    			}else{
        			GeneExp <- merge(GeneExp, bindtemp, by="GENE", sort=F)
			}	
		}
	}

	save(GeneExp, file = paste(path, "GeneExp.rda", sep="/"))

	if (outputList){
		write.csv(GeneExp, file = paste(path, "GeneExp.csv", sep="/"), row.names=F)
	}
}