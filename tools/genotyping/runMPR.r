#!/opt/blc/genome/biosoft/R/bin/Rscript

self='./runMPR.r';
argv = commandArgs(T);

if (is.null(argv) | length(argv)<2) {
  cat("Error: No options. Run",self,"<input> <output_prefix>.{pGT,rgt,bin,png,block}\n")
  q(status=1)
}

sink(paste(argv[2],'.log',sep=''))

source('./mpr.R');

InDat=read.table(argv[1],T,na.strings='-',comment.char = '',as.is=TRUE)
#InDat=read.table('psnp/Chr01.add_ref',T,na.strings='-',comment.char = '',as.is=TRUE)
Ref=data.frame(InDat$Ref,row.names=InDat$Pos,check.names=F)
SNP=as.matrix(data.frame(InDat[-(1:3)],row.names=InDat$Pos))
#InDat='';

#myBaseData <- SNP[sample(120,50),]
myBaseData <- SNP
snpT=base2Allele(myBaseData)
Tmp=sort(as.numeric(rownames(snpT)))
myBaseData <- myBaseData[as.character(Tmp),]

allele.MPR <- localMPR(baseData=myBaseData,maxIterate=50,maxNStep=5,showDetail=TRUE,verbose=TRUE)
warnings()

#snpSet <- sort(as.numeric(rownames(na.omit(allele.MPR))))
#allele.MPR <- allele.MPR[match(snpSet,rownames(allele.MPR)),]

write.table(allele.MPR,paste(argv[2],'.r0pGT',sep=''),col.names=TRUE,quote=F,sep='\t')

system.time(all.res <- globalMPRRefine(myBaseData,alleleA=allele.MPR[,1],groupSort=TRUE,
             numPerm=20,numTry=3,
             numBaseStep=50,numBaseCandidateStep=100,numKnownStep=30,
             numKnownCandidateStep=50,useMedianToFindKnown=TRUE,maxIterate=150,
             maxNStep=3,scoreMin=0.8,saveMidData=TRUE,verbose=TRUE))
warnings()

write.table(allele.MPR,paste(argv[2],'.r1pGT',sep=''),col.names=TRUE,quote=F,sep='\t')

perm <- 10
res <- all.res$midData[[perm]]
table(geno.res <- genotypeCallsBayes(res$call,errorRate=1e-11,eps=1e-10,maxIterate=100,verbose=FALSE)$type)

allele.MPR <- res$allele
allele.MPR[geno.res==2|rowMin(res$call)>=perm,] <- NA
allele.MPR[geno.res==3,] <- allele.MPR[geno.res==3,c(2,1)]

table(is.na(alleleA <- allele.MPR[,1]))
#table(is.na(alleleB <- allele.MPR[,2]))
idsA <- match(names(alleleA),rownames(Ref))
#idsB <- match(names(alleleB),rownames(Ref))
table(Ref[idsA,1]==alleleA)
#table(Ref[idsB,1]==alleleB)

snpSet <- sort(as.numeric(rownames(na.omit(allele.MPR))))
#alleleMPR=allele.MPR[match(snpSet,rownames(allele.MPR)),]

#argv=c('iChr01','oChr01')

write.table(allele.MPR,paste(argv[2],'.pGT',sep=''),col.names=TRUE,quote=F,sep='\t')

ids <- match(snpSet,rownames(myBaseData))
table(is.na(geno.data <- base2Geno(myBaseData[ids,],allele.MPR[ids,])))
write.table(geno.data,paste(argv[2],'.rgt',sep=''),col.names=TRUE,quote=F,sep='\t')

theBin=genoToBin(geno.data,as.numeric(rownames(geno.data)),corrected=F,heterozygote=TRUE,minBinsize=250000,geno.probability=c(0.49951171875, 0.49951171875,0.0009765625))
warnings()

write.table(theBin$bin,paste(argv[2],'.bin',sep=''),col.names=TRUE,quote=F,sep='\t')
write.table(theBin$border,paste(argv[2],'.border',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
write.table(theBin$block$'1'$block,paste(argv[2],'.block',sep=''),col.names=F,quote=F,sep='\t')

geno.bin=theBin$bin

geno.colors <- geno.bin;geno.colors[is.na(geno.colors)] <- rgb(0,0,0)
geno.colors[geno.colors==0]=rgb(0,0.7,0)
geno.colors[geno.colors==1]=rgb(0,0,0.7)
geno.colors[geno.colors==0.5]=rgb(0.6,0,0)
rils=matrix(theBin$border,nrow=nrow(geno.bin),ncol=ncol(geno.bin))	# 1:nrow(geno.bin) or theBin$border
poses=matrix(rep(1:ncol(geno.bin),each=nrow(geno.bin)),nrow=nrow(geno.bin),ncol=ncol(geno.bin))
png(paste(argv[2],'.png',sep=''), 1200, 4800, bg='white')
plot(poses, rils, col=geno.colors,pch=20, ylab='',xlab="RIL index",mar=c(1.1,2.1,2.1,1.1))

for(j in 1:ncol(geno.colors)) {
	for(i in 1:(nrow(geno.colors)-1)) {
		if(geno.colors[i,j]==geno.colors[i+1,j])
			segments(j,rils[i,j],j,rils[i+1,j],col=geno.colors[i,j],lwd=2)
	}
}

for(i in 1:nrow(rils)) {
	segments(1,rils[i,1],ncol(rils),rils[i,ncol(rils)],lwd=1,col=rgb(0,0,0.2,0.4))
}

dev.off()
