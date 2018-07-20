#!/usr/bin/env Rscript

# install.packages('devtools')
# devtools::install_github('alyssafrazee/RSkittleBrewer')
# source("https://bioconductor.org/biocLite.R")
# biocLite("ballgown")

library(plyr)

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

## Read phenotype sample data
#pheno_data <- read.csv('sample.csv')
pheno_data <- data.frame(row.names = c(1,2,3),ids=c('expA','expB','wild'),kind=c('Exp','Exp','Wild'))

bg_chrX <- ballgown(dataDir = "out", samplePattern="", pData=pheno_data)
bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)
#results_transcripts <-  stattest(bg_chrX_filt, feature='transcript', covariate="kind", getFC=TRUE, meas='FPKM')
#results_genes <-  stattest(bg_chrX_filt, feature='gene', covariate='kind', getFC=TRUE, meas='FPKM')

## Add gene name
#results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

## Sort results from smallest p-value
#results_transcripts <- arrange(results_transcripts, pval)
#results_genes <-  arrange(results_genes, pval)

## Write results to CSV
#write.csv(results_transcripts, "new_transcripts_results.csv", row.names=FALSE)
#write.csv(results_genes, "new_genes_results.csv", row.names=FALSE)

bg <- bg_chrX_filt
#bg <- bg_chrX
# whole_tx_table
wtx = texpr(bg, 'all')
# write.csv(wtx,'wtx.csv', row.names=FALSE)

res <- cbind(wtx,'fpkmAtW'=log2(wtx$FPKM.expA/wtx$FPKM.wild),'fpkmBtW'=log2(wtx$FPKM.expB/wtx$FPKM.wild),'fpkmDtA'=2*(wtx$FPKM.expA-wtx$FPKM.expB)/(wtx$FPKM.expA+wtx$FPKM.expB))
#x <- colnames(res)
#x[length(x)-2]<-'fpkmAtW'
#x[length(x)-1]<-'fpkmBtW'
#x[length(x)]<-'fpkmDtA'
#colnames(res) <- x
#head(res)
#res <- na.omit(res)

t <- cor(res$FPKM.expA,res$FPKM.expB,use='pairwise.complete.obs',method='pearson')
cat("Corr:",t,"\n")

getype <- function(Vx,Vy) {
	theType = 'Z'
	if (is.na(Vx) || is.na(Vy) || is.na(Vx-Vy)) {
		theType = 'N'
	} else {
		if (abs(Vx)<1 & abs(Vy)<1) {
			theType = 'O'
		} else if (abs(Vx-Vy)<=1 & (abs(Vx)>=1 | abs(Vy)>=1)) {
			theType ='A'
		} else if ((abs(Vx)<1 & abs(Vy)>1) | (abs(Vy)<1 & abs(Vx)>1)) {
			theType ='B'
		} else if (abs(Vx)>1 | abs(Vy)>1) {
			if (Vx*Vy<0) {
				theType ='C'
			} else if (abs(Vx)>abs(Vy)) {
				theType ='D'
			} else if (abs(Vx)<abs(Vy)) {
				theType ='E'
			} else {
				theType ='X'
			}
		}
	}
	theType
}
#apply(head(res)[,c('fpkmAtW','fpkmBtW')], 1, function(y) getype(y['fpkmAtW'],y['fpkmBtW']))
#df=head(res)
#cbind(df,'Type' = mapply(getype, df$fpkmAtW, df$fpkmBtW) )
res <- cbind(res,'Type' = mapply(getype, res$fpkmAtW, res$fpkmBtW) )
#head(res[res$Type=='Z',])

#write.csv(res,'res.csv', row.names=FALSE)

Vx <- res$fpkmAtW
Vy <- res$fpkmBtW
ex_input_p_mvc <- rbind('Vx'=res$fpkmAtW,'Vy'=res$fpkmBtW)

mvc_cross<-ex_input_p_mvc[,((abs(Vx)<1 & abs(Vy)>1) | (abs(Vy)<1 & abs(Vx)>1))]
mvc_trend<-ex_input_p_mvc[,(abs(Vx-Vy)<1 & (abs(Vx)>1 | abs(Vy)>1))]
mvc_gray<-ex_input_p_mvc[,(abs(Vx)>1 | abs(Vy)>1)]

pdf('res.pdf',title='NSFC Report Fig.1',pointsize=16)

plot(mvc_gray[1,],mvc_gray[2,],xlim=c(-20,20),ylim=c(-20,20),pch='.',cex=3,
xlab='Relative expression level in mutA SW-13-2-B7A (log2ratio)',
ylab='Relative expression level in mutB SW-13-2-D3A (log2ratio)',
main='Deregulated mRNAs',col='darkgray')

points(mvc_trend[1,],mvc_trend[2,],pch='.',cex=3,col='blue')
points(mvc_cross[1,],mvc_cross[2,],pch='.',cex=3,col='red')

abline(h=0,col='black',lty=4)
abline(v=0,col='black',lty=4)

Ldis <- 17
Sdis <- 12
Tdis <- 5
text(c(-Ldis,Ldis),c(-Ldis,Ldis),"A",cex=2)
text(c(0,0),c(-Ldis,Ldis),"B",cex=2)
text(c(Sdis,-Sdis),c(-Sdis,Sdis),"C",cex=2)
text(c(Tdis,-Tdis),c(Sdis,-Sdis),"E",cex=2)
text(c(Sdis,-Sdis),c(Tdis,-Tdis),"D",cex=2)

dev.off()

#############3
#slotNames(bg_chrX_filt)
#head(bg_chrX_filt@structure)

library('rtracklayer')
gtf<-gzfile('gencode.v28.primary_assembly.annotation.gtf.gz','rt')
rawgtf <- import(gtf)
mygtf <- rawgtf[rawgtf$type=='gene']

stra <- revalue(res$strand, c('.' = "*"))
resGR <- GRanges(Rle(res$chr),IRanges(res$start,res$end,names=res$t_id),Rle(strand(stra)),
                 FPKM.expA=res$FPKM.expA, FPKM.expB=res$FPKM.expB, FPKM.wild=res$FPKM.wild)

overlapGenes <- findOverlaps(resGR, mygtf)
overlapGenes.df <- as.data.frame(overlapGenes)
ohit=mygtf[overlapGenes.df$subjectHits]
#head(ohit$gene_id)
ores <- cbind(res[overlapGenes.df$queryHits,],'GeneID'=ohit$gene_id,'GeneType'=ohit$gene_type,'GeneName'=ohit$gene_name)
ores <- ores[order(ores$Type),]
#nearestGenes <- distanceToNearest(resGR, mygtf)
#nearestGenes.df <- as.data.frame(overlapGenes)
#nhit=mygtf[nearestGenes.df$subjectHits]
#nres <- cbind(res[nearestGenes.df$queryHits,],nhit$gene_id,nhit$gene_type,nhit$gene_name)

#x <- colnames(ores)
#x[length(x)-2]<-'GeneID'
#x[length(x)-1]<-'GeneType'
#x[length(x)]<-'GeneName'
#colnames(ores) <- x
write.csv(ores,'ores.csv', row.names=FALSE)
