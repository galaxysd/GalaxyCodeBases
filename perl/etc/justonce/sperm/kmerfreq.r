#!/usr/bin/env r
if (is.null(argv) | length(argv)<3) {
  cat("Usage: kmerfreq.r maxcnt hist_file coverage [scaledmin] [scaledmax]\n")
  q()
}
testXrange <- as.numeric(argv[1])
infile <- argv[2]
avgcvg <- as.numeric(argv[3])
scaledmin <- as.numeric(argv[4])
scaledmax <- as.numeric(argv[5])

if(!exists("D__myfit"))
	source("myfit.R")

#inid <- unlist(strsplit(infile,'.',TRUE))[1]
inid <- strsplit(infile,'.',TRUE)[[1]][1]
TheIDRelshp <- c('Donor','MALBAC Sperm01','MALBAC Sperm02','MALBAC Sperm03','MDA Sperm23','MDA Sperm24','MDA Sperm28','Donor')
names(TheIDRelshp) <- c('blood','S1','S2','S3','S23','S24','S28','mlbacDonor')

maxXrange <- 25
"
testXrange <- 15
infile <- 'blood.fq.k25.gz.hist.fixed'
avgcvg <- 3.27885
"
library('vcd')

aa=read.table(infile,skip=8)
thefreq <- aa[,1]
cnt <- as.numeric(aa[,2])
theratio <- aa[,3]

thedata <- cbind(thefreq,cnt)
flag <- thedata[,1] <= testXrange
#usedata <- cbind(cnt[flag],thefreq[flag])
usedata <- cbind(theratio[flag],thefreq[flag])

scaled <- scale2pois(usedata,avgcvg,cntmin=scaledmin,cntmax=scaledmax)
scaledata <- scaled$ret
scaleres <- myfit(scaledata, par = list(lambda = avgcvg))
print('===scaledmyfit===')
#print(scaleres)
scalesum <- summary(scaleres)

scalechi <- chisq.test(cbind(scaleres$observed,scaleres$fitted))
print(scalechi)

fitres <- goodfit(usedata,type='poisson', par = list(lambda = avgcvg))
#fitres <- goodfit(usedata,type='poisson', method = 'MinChisq')
#sumres <- summary(fitres)

selfres <- myfit(usedata, par = list(lambda = avgcvg))
print('===myfit===')
#print(selfres)
selfsum <- summary(selfres)
#print(selfsum)

newres<-list(observed = fitres$observed[c(-1,-2)], count = fitres$count[c(-1,-2)], fitted = fitres$fitted[c(-1,-2)],
		type = "poisson", method = fitres$method, df = fitres$df - 2, par = fitres$par)
class(newres) <- "goodfit"
print('===goodfit===')
#print(newres)
sumres <- summary(newres)
#print(sumres)

dbg1 <- function(x) {
	arg <- deparse(substitute(x))
	RVAL <- cbind(x$observed,x$fitted,(x$observed-x$fitted)/x$fitted)
	#RVAL <- round(RVAL,digits = 6)
	colnames(RVAL) <- c(paste0(arg,'.ob'),paste0(arg,'.fi'),'(O-E)/E')
	rownames(RVAL) <- 1+ 1:nrow(RVAL)
	RVAL
}
print(cbind(dbg1(scaleres),dbg1(newres),dbg1(selfres)))

tiff(filename = paste0(infile,".tiff"), compression="lzw", width = 683, height = 683)#, units = "px", pointsize = 52)
par(mar=c(5, 6, 2, 2),ps=20,family='sans')

yy <- theratio
xx <- thefreq

mean2=avgcvg
xxx=seq(0,max(xx))
zz=dpois(xxx,mean2);
maxY=max(zz,yy,zz*selfres$zoom)
maxY<-0.32
#lines(xx,zz,xlim=c(0,60),type='l');

#barplot(yy,xlim=c(0,20),ylim=c(0,0.212),xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
plot(xx,yy,type='h',lwd=19,xlim=c(1,maxXrange),ylim=c(0,maxY),
	main = paste0("Pearson goodness of fit, p=",signif(scalesum[1,3],digits = 6)),
	xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
# dpois(3,3.62)=0.2117524
legend("topright",pch=-1,lty=1,col=c('navy','green3','red'), x.intersp = 1, y.intersp = 2,cex=1,lwd=c(11,7,8),
	legend= c( TheIDRelshp[inid],paste0("lambda=",signif(avgcvg,digits=6)),
		paste0("scaledby=",signif(scaled$zoom, digits=6)) )
)
lines(scaleres$count,scaleres$observed,xlim=c(0,60),type='h',col=rgb(1,1,0,1),lwd=11)

axis(1, at = seq(0, maxXrange, by = 5),lwd=3,cex=1)
lines(xxx,zz,xlim=c(0,maxXrange),type='l',col=rgb(0,0.804,0,0.8),lwd=7)
#lines(xx,yy*scaled$zoom,xlim=c(0,60),type='l',col=rgb(1,0,0,0.7),lwd=3)
#lines(scaleres$count,scaleres$observed,xlim=c(0,60),type='l',col=rgb(1,0,0,0.9),lwd=5)

lines(xx,yy*scaled$zoom,xlim=c(0,60),type='h',col=rgb(1,0,0,1),lwd=7)
dev.off()

print(selfsum)
outline <- paste0(TheIDRelshp[inid],':[',scaled$cntmin,',',scaled$cntmax,"] cvg=",avgcvg," scaleR=",scaled$zoom)
print(outline)
write(outline,paste0(infile,".txt"),append=FALSE)

outline <- paste0('X^2= ',scalesum[1,1],', df= ',scalesum[1,2],', p= ',scalesum[1,3])
print(outline)
write(outline,paste0(infile,".txt"),append=TRUE)

print(scalesum)

"
#./kmerfreq.r mlbacDonor.k25.lz4.hist 6.85
#./kmerfreq.r mdaS23.k25.lz4.hist 4.34
#./kmerfreq.r S01.k25.lz4.hist 3.62

#                         X^2 df  P(> X^2)
#Pearson           1.56832025 18 0.9999998
#Likelihood Ratio -0.06964308 18 1.0000000

#                          X^2 df P(> X^2)
#Pearson          7285.4140548 18        0
#Likelihood Ratio    0.3987449 18        1

#                          X^2 df  P(> X^2)
#Pearson          2.511278e+05 18 0.0000000
#Likelihood Ratio 2.540676e+00 18 0.9999924

./kmerfreq.r 15 mlbacDonor.k25.lz4.hist.fixed 5.67811
./kmerfreq.r 15 S1.fq.k25.gz.hist.fixed 3.29571
./kmerfreq.r 15 S2.fq.k25.gz.hist.fixed 3.28174
./kmerfreq.r 15 S3.fq.k25.gz.hist.fixed 2.94459
./kmerfreq.r 15 S23.fq.k25.gz.hist.fixed 3.26044
./kmerfreq.r 15 S24.fq.k25.gz.hist.fixed 3.27885
./kmerfreq.r 15 S28.fq.k25.gz.hist.fixed 3.30884
./kmerfreq.r 15 blood.fq.k25.gz.hist.fixed 3.27885

"
