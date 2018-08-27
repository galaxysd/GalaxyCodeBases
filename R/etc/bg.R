#!/usr/bin/env Rscript
# run this in the output directory for rnaseq_pipeline.sh
# passing the pheno data csv file as the only argument 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
# assume no output directory argument was given to rnaseq_pipeline.sh
  pheno_data_file <- paste0(getwd(),"/sample.csv")
} else {
  pheno_data_file <- args[1]
}

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

## Read phenotype sample data
pheno_data <- read.csv(pheno_data_file)

## Read in expression data
out <- ballgown(dataDir = "out", samplePattern= "", pData=pheno_data)

## Filter low abundance genes
out_filt <- subset(out, "rowVars(texpr(out)) > 1", genomesubset=TRUE)

## Divided by group
sampleL <- subset(out_filt,"tissue %in% 'L'", genomesubset=FALSE)
sampleB <- subset(out_filt,"tissue %in% 'B'", genomesubset=FALSE)
sampleE <- subset(out_filt,"tissue %in% 'E'", genomesubset=FALSE)

## Divided by content
LA <- subset(sampleL,"gko %in% 'WT'", genomesubset=FALSE)
LB <- subset(sampleL,"exp %in% 'CD'", genomesubset=FALSE)
LC <- subset(sampleL,"tko %in% c('CD-WT','HFD-KO')", genomesubset=FALSE)
LD <- subset(sampleL,"exp %in% 'HFD'", genomesubset=FALSE)
BA <- subset(sampleB,"gko %in% 'WT'", genomesubset=FALSE)
BB <- subset(sampleB,"exp %in% 'CD'", genomesubset=FALSE)
BC <- subset(sampleB,"tko %in% c('CD-WT','HFD-KO')", genomesubset=FALSE)
BD <- subset(sampleB,"exp %in% 'HFD'", genomesubset=FALSE)
EA <- subset(sampleE,"gko %in% 'WT'", genomesubset=FALSE)
EB <- subset(sampleE,"exp %in% 'CD'", genomesubset=FALSE)
EC <- subset(sampleE,"tko %in% c('CD-WT','HFD-KO')", genomesubset=FALSE)
ED <- subset(sampleE,"exp %in% 'HFD'", genomesubset=FALSE)
LA_data <- subset(pheno_data,tissue=="L"&gko=="WT",genomesubset=FALSE)

## DE by transcript
LAresults_transcripts <-  stattest(LA, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
LBresults_transcripts <-  stattest(LB, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
LCresults_transcripts <-  stattest(LC, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
LDresults_transcripts <-  stattest(LD, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
BAresults_transcripts <-  stattest(BA, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
BBresults_transcripts <-  stattest(BB, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
BCresults_transcripts <-  stattest(BC, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
BDresults_transcripts <-  stattest(BD, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
EAresults_transcripts <-  stattest(EA, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
EBresults_transcripts <-  stattest(EB, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
ECresults_transcripts <-  stattest(EC, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')
EDresults_transcripts <-  stattest(ED, feature='transcript', covariate="tko", getFC=TRUE, meas='FPKM')

## DE by gene
LAresults_genes <- stattest(LA, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
LBresults_genes <- stattest(LB, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
LCresults_genes <- stattest(LC, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
LDresults_genes <- stattest(LD, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
BAresults_genes <- stattest(BA, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
BBresults_genes <- stattest(BB, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
BCresults_genes <- stattest(BC, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
BDresults_genes <- stattest(BD, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
EAresults_genes <- stattest(EA, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
EBresults_genes <- stattest(EB, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
ECresults_genes <- stattest(EC, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')
EDresults_genes <- stattest(ED, feature='gene', covariate='tko',getFC=TRUE, meas='FPKM')

## Add gene name
LAresults_transcripts <- data.frame(chr=LA@expr$trans$chr,strand=LA@expr$trans$strand,start=LA@expr$trans$start,end=LA@expr$trans$end,geneNames=ballgown::geneNames(LA),geneIDs=ballgown::geneIDs(LA), LAresults_transcripts)
LBresults_transcripts <- data.frame(chr=LB@expr$trans$chr,strand=LB@expr$trans$strand,start=LB@expr$trans$start,end=LB@expr$trans$end,geneNames=ballgown::geneNames(LB),geneIDs=ballgown::geneIDs(LB), LBresults_transcripts)
LCresults_transcripts <- data.frame(chr=LC@expr$trans$chr,strand=LC@expr$trans$strand,start=LC@expr$trans$start,end=LC@expr$trans$end,geneNames=ballgown::geneNames(LC),geneIDs=ballgown::geneIDs(LC), LCresults_transcripts)
LDresults_transcripts <- data.frame(chr=LD@expr$trans$chr,strand=LD@expr$trans$strand,start=LD@expr$trans$start,end=LD@expr$trans$end,geneNames=ballgown::geneNames(LD),geneIDs=ballgown::geneIDs(LD), LDresults_transcripts)
BAresults_transcripts <- data.frame(chr=BA@expr$trans$chr,strand=BA@expr$trans$strand,start=BA@expr$trans$start,end=BA@expr$trans$end,geneNames=ballgown::geneNames(BA),geneIDs=ballgown::geneIDs(BA), BAresults_transcripts)
BBresults_transcripts <- data.frame(chr=BB@expr$trans$chr,strand=BB@expr$trans$strand,start=BB@expr$trans$start,end=BB@expr$trans$end,geneNames=ballgown::geneNames(BB),geneIDs=ballgown::geneIDs(BB), BBresults_transcripts)
BCresults_transcripts <- data.frame(chr=BC@expr$trans$chr,strand=BC@expr$trans$strand,start=BC@expr$trans$start,end=BC@expr$trans$end,geneNames=ballgown::geneNames(BC),geneIDs=ballgown::geneIDs(BC), BCresults_transcripts)
BDresults_transcripts <- data.frame(chr=BD@expr$trans$chr,strand=BD@expr$trans$strand,start=BD@expr$trans$start,end=BD@expr$trans$end,geneNames=ballgown::geneNames(BD),geneIDs=ballgown::geneIDs(BD), BDresults_transcripts)
EAresults_transcripts <- data.frame(chr=EA@expr$trans$chr,strand=EA@expr$trans$strand,start=EA@expr$trans$start,end=EA@expr$trans$end,geneNames=ballgown::geneNames(EA),geneIDs=ballgown::geneIDs(EA), EAresults_transcripts)
EBresults_transcripts <- data.frame(chr=EB@expr$trans$chr,strand=EB@expr$trans$strand,start=EB@expr$trans$start,end=EB@expr$trans$end,geneNames=ballgown::geneNames(EB),geneIDs=ballgown::geneIDs(EB), EBresults_transcripts)
ECresults_transcripts <- data.frame(chr=EC@expr$trans$chr,strand=EC@expr$trans$strand,start=EC@expr$trans$start,end=EC@expr$trans$end,geneNames=ballgown::geneNames(EC),geneIDs=ballgown::geneIDs(EC), ECresults_transcripts)
EDresults_transcripts <- data.frame(chr=ED@expr$trans$chr,strand=ED@expr$trans$strand,start=ED@expr$trans$start,end=ED@expr$trans$end,geneNames=ballgown::geneNames(ED),geneIDs=ballgown::geneIDs(ED), EDresults_transcripts)

## Sort results from smallest p-value
LAresults_transcripts <- arrange(LAresults_transcripts, pval)
LAresults_genes <-  arrange(LAresults_genes, pval)
LBresults_transcripts <- arrange(LBresults_transcripts, pval)
LBresults_genes <-  arrange(LBresults_genes, pval)
LCresults_transcripts <- arrange(LCresults_transcripts, pval)
LCresults_genes <-  arrange(LCresults_genes, pval)
LDresults_transcripts <- arrange(LDresults_transcripts, pval)
LDresults_genes <-  arrange(LDresults_genes, pval)
BAresults_transcripts <- arrange(BAresults_transcripts, pval)
BAresults_genes <-  arrange(BAresults_genes, pval)
BBresults_transcripts <- arrange(BBresults_transcripts, pval)
BBresults_genes <-  arrange(BBresults_genes, pval)
BCresults_transcripts <- arrange(BCresults_transcripts, pval)
BCresults_genes <-  arrange(BCresults_genes, pval)
BDresults_transcripts <- arrange(BDresults_transcripts, pval)
BDresults_genes <-  arrange(BDresults_genes, pval)
EAresults_transcripts <- arrange(EAresults_transcripts, pval)
EAresults_genes <-  arrange(EAresults_genes, pval)
EBresults_transcripts <- arrange(EBresults_transcripts, pval)
EBresults_genes <-  arrange(EBresults_genes, pval)
ECresults_transcripts <- arrange(ECresults_transcripts, pval)
ECresults_genes <-  arrange(ECresults_genes, pval)
EDresults_transcripts <- arrange(EDresults_transcripts, pval)
EDresults_genes <-  arrange(EDresults_genes, pval)

## Write results to CSV
write.csv(LAresults_transcripts, "LA_transcripts_results.csv", row.names=FALSE)
write.csv(LAresults_genes, "LA_genes_results.csv", row.names=FALSE)
write.csv(LBresults_transcripts, "LB_transcripts_results.csv", row.names=FALSE)
write.csv(LBresults_genes, "LB_genes_results.csv", row.names=FALSE)
write.csv(LCresults_transcripts, "LC_transcripts_results.csv", row.names=FALSE)
write.csv(LCresults_genes, "LC_genes_results.csv", row.names=FALSE)
write.csv(LDresults_transcripts, "LD_transcripts_results.csv", row.names=FALSE)
write.csv(LDresults_genes, "LD_genes_results.csv", row.names=FALSE)
write.csv(BAresults_transcripts, "BA_transcripts_results.csv", row.names=FALSE)
write.csv(BAresults_genes, "BA_genes_results.csv", row.names=FALSE)
write.csv(BBresults_transcripts, "BB_transcripts_results.csv", row.names=FALSE)
write.csv(BBresults_genes, "BB_genes_results.csv", row.names=FALSE)
write.csv(BCresults_transcripts, "BC_transcripts_results.csv", row.names=FALSE)
write.csv(BCresults_genes, "BC_genes_results.csv", row.names=FALSE)
write.csv(BDresults_transcripts, "BD_transcripts_results.csv", row.names=FALSE)
write.csv(BDresults_genes, "BD_genes_results.csv", row.names=FALSE)
write.csv(EAresults_transcripts, "EA_transcripts_results.csv", row.names=FALSE)
write.csv(EAresults_genes, "EA_genes_results.csv", row.names=FALSE)
write.csv(EBresults_transcripts, "EB_transcripts_results.csv", row.names=FALSE)
write.csv(EBresults_genes, "EB_genes_results.csv", row.names=FALSE)
write.csv(ECresults_transcripts, "EC_transcripts_results.csv", row.names=FALSE)
write.csv(ECresults_genes, "EC_genes_results.csv", row.names=FALSE)
write.csv(EDresults_transcripts, "ED_transcripts_results.csv", row.names=FALSE)
write.csv(EDresults_genes, "ED_genes_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
subset(LAresults_transcripts, LAresults_transcripts$qval <=0.05)
subset(LAresults_genes, LAresults_genes$qval <=0.05)

## Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
#fpkm <- texpr(bg_chrX, meas='FPKM')
#fpkm <- log2(fpkm +1)
#boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')

## Plot individual transcripts
#ballgown::transcriptNames(bg_chrX)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_chrX)[12], ' : ',ballgown::transcriptNames(bg_chrX)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
#plotMeans(ballgown::geneIDs(bg_chrX)[203], bg_chrX_filt, groupvar="sex", legend=FALSE)

##print gene abundance distribution in your screen
pdf(file="Gene_abundance_distribution.pdf",width=15,height=9)
par(mai=c(1.8,1,1,1))
tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
fpkm <- texpr(out_filt, meas='FPKM')
fpkm <- log2(fpkm +1)
boxplot(fpkm, col=as.numeric(pheno_data$tko), las=2,ylab='log2(FPKM+1)')
dev.off()

##print individual transcripts
pdf(file="Individual_transcripts.pdf")
ballgown::transcriptNames(out_filt)[15]
plot(fpkm[15,] ~ pheno_data$gko, border=c(1,2),main=paste(ballgown::geneNames(out_filt)[15], ' : ',ballgown::transcriptNames(out_filt)[15]),pch=19, xlab="Gko", ylab='log2(FPKM+1)')
points(fpkm[15,] ~ jitter(as.numeric(pheno_data$gko)), col=as.numeric(pheno_data$gko))
dev.off()

##print Plot gene of transcript
pdf(file="Plot_gene_of_transcript.pdf")
plotTranscripts(ballgown::geneIDs(out_filt)[8219], out_filt,main=c('Gene XIST in sample HFD-L14'), sample=c('HFD-L14'))
dev.off()

##print the compare of average expression
pdf(file="Comparison_of_average_expression.pdf")
plotTranscripts(ballgown::geneIDs(out_filt)[8219], out_filt,main=c('Gene XIST in sample HFD-L14'), sample=c('HFD-L14'))
plotMeans('MSTRG.7196', out_filt,groupvar="gko",legend=FALSE)
dev.off()

#print the first picture in sequence
fpkm2 <- texpr(out_filt, meas='FPKM')
tbl_fpkm2<-tbl_df(fpkm2)
colnames(tbl_fpkm2) <- c("CD-B1","CD-B2","CD-B3","CD-B4","CD-B5","CD-E1","CD-E2","CD-E4","CD-E5","CD-KO-B10","CD-KO-B6","CD-KO-B7","CD-KO-B8","CD-KO-E10","CD-KO-E6","CD-KO-E7","CD-KO-E8","CD-KO-E9","CD-KO-L10","CD-KO-L6","CD-KO-L7","CD-KO-L8","CD-KO-L9","CD-L1","CD-L2","CD-L3","CD-L4","CD-L5","HFD-B11","HFD-B12","HFD-B13","HFD-B14","HFD-B15","HFD-E11","HFD-E12","HFD-E13","HFD-E14","HFD-E15","HFD-KO-B16","HFD-KO-B17","HFD-KO-B18","HFD-KO-B19","HFD-KO-E16","HFD-KO-E17","HFD-KO-E18","HFD-KO-E19","HFD-KO-E20","HFD-KO-L16","HFD-KO-L17","HFD-KO-L18","HFD-KO-L19","HFD-KO-L20","HFD-L11","HFD-L12","HFD-L13","HFD-L14","HFD-L15")
newdf <- select(tbl_fpkm2,"CD-L1","CD-L2","CD-L3","CD-L4","CD-L5","CD-KO-L6","CD-KO-L7","CD-KO-L8","CD-KO-L9","CD-KO-L10","HFD-L11","HFD-L12","HFD-L13","HFD-L14","HFD-L15","HFD-KO-L16","HFD-KO-L17","HFD-KO-L18","HFD-KO-L19","HFD-KO-L20","CD-B1","CD-B2","CD-B3","CD-B4","CD-B5","CD-KO-B6","CD-KO-B7","CD-KO-B8","CD-KO-B10","HFD-B11","HFD-B12","HFD-B13","HFD-B14","HFD-B15","HFD-KO-B16","HFD-KO-B17","HFD-KO-B18","HFD-KO-B19","CD-E1","CD-E2","CD-E4","CD-E5","CD-KO-E6","CD-KO-E7","CD-KO-E8","CD-KO-E9","CD-KO-E10","HFD-E11","HFD-E12","HFD-E13","HFD-E14","HFD-E15","HFD-KO-E16","HFD-KO-E17","HFD-KO-E18","HFD-KO-E19","HFD-KO-E20")
pdf(file="Gene_abundance_distribution_insequnce.pdf",width=15,height=9)
par(mai=c(1.8,1,1,1))
tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
newdf <- log2(newdf +1)
boxplot(newdf, col=as.numeric(pheno_data$tko), las=2,ylab='log2(FPKM+1)')
dev.off()

