source('./mpr.R');

InDat=read.table('psnp/Chr01.add_ref',T,na.strings='-',comment.char = '',as.is=T)
Ref=data.frame(InDat$Ref,row.names=InDat$Pos,check.names=F)
SNP=as.matrix(data.frame(InDat[-(1:3)],row.names=InDat$Pos))
#InDat='';

myBaseData <- SNP[sample(120,50),]


allele.MPR <- localMPR(baseData=myBaseData,maxIterate=50,maxNStep=5,showDetail=TRUE)

snpSet <- sort(as.numeric(rownames(na.omit(allele.MPR))))
allele.MPR <- allele.MPR[match(snpSet,rownames(allele.MPR)),]

system.time(all.res <- globalMPRRefine(myBaseData,alleleA=na.omit(
             allele.MPR[,1]),groupSort=TRUE,numPerm=20,numTry=3,
             numBaseStep=50,numBaseCandidateStep=100,numKnownStep=30,
             numKnownCandidateStep=50,useMedianToFindKnown=TRUE,maxIterate=150,
             maxNStep=3,scoreMin=0.8,saveMidData=TRUE,verbose=TRUE))

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
ids <- match(snpSet,rownames(myBaseData))
table(is.na(geno.data <- base2Geno(myBaseData[ids,],allele.MPR[ids,])))

