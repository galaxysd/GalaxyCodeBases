#!/usr/bin/env littler

if (length(argv)<2) {
    cat("Usage ./dodseq.r <prefix.featureCounts.txt> <Nctl,Nexp> [prefix]\n")
    quit()
}
inFile <- argv[1]
ExpCnt <- as.numeric(unlist(strsplit(argv[2],',',fixed=T)))
outP <- argv[3]
if (length(argv)==2) {
    outP <- unlist(strsplit(inFile,'.',fixed=T))[1]
}
cat(sep='',"[!]From [",inFile,"],Design [",ExpCnt[1],':',ExpCnt[2],"],Prefix [",outP,"].\n")

countdata <- read.table(inFile, header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("^alnSTAT\\.", "", colnames(countdata))
colnames(countdata) <- gsub("Aligned\\.sortedByCoord\\.out\\.[sb]am$", "", colnames(countdata))
countdata <- as.matrix(countdata)

tFileDim <- dim(countdata)
tExpSize <- ExpCnt[1]+ExpCnt[2]
if ( tFileDim[2] != tExpSize ) {
    cat("[x]Design size",tExpSize,"not matching with file column size",tFileDim[2],"!\n[!]File columns:[")
    cat(sep=',',colnames(countdata))
    cat("].\n")
    quit()
}

cat("[!]File Data preview:\n")
print(head(countdata))

(condition <- factor(c(rep("ctl", ExpCnt[1]), rep("exp", ExpCnt[2]))))
#print(condition)

library(DESeq2)
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
