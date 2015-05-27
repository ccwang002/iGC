hg19DBALL <- read.csv("C:/Users/Daphne/Documents/Project/hg19-RefSeq.csv", header = TRUE)

index <- grep("NM", hg19DBALL[,7])
hg19DBNM = hg19DBALL[index,1:7]

save(hg19DBNM,file="C:/Users/Daphne/Documents/Project/hg19DBNM.rda")
write.csv(hg19DBNM, "hg19DBNM.csv", row.names=F)

WG <- unique(hg19DBNM[,6])
save(WG, file="C:/Users/Daphne/Documents/Project/wholegenes.rda")
