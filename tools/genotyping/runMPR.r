source('./mpr.R');

InDat=read.table('psnp/Chr01.add_ref',T,na.strings='-',comment.char = '',as.is=T)
Ref=data.frame(InDat$Ref,row.names=InDat$Pos,check.names=F)
SNP=data.frame(InDat[-(1:3)],row.names=InDat$Pos)
#InDat='';


if (length(b[b %in% "W"])) {cat("Done.\n")}


allele.MPR <- localMPR(baseData=SNP,maxIterate=50,maxNStep=5,showDetail=TRUE)








myBaseData <- snpData[sample(200,50),]

## MPR inference with maximum step size of 5
allele.MPR <- localMPR(baseData=myBaseData,maxIterate=50,maxNStep=5,showDetail=TRUE)
