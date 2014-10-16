#!/usr/bin/env r
if (is.null(argv) | length(argv)<2) {
  cat("Usage: kmerfreq.r blood_hist hist_file\n")
  q()
}
bloodfile <- argv[1]
infile <- argv[2]

if(!exists("D__myfit"))
	source("myfit.R")

maxXrange <- 30
testXrange <- 15
#infile <- 'mlbacDonor.k25.lz4.hist'
#avgcvg <- 6.85

#library('vcd')

blood=read.table(bloodfile,skip=8)
#mda=read.table('mdaS23.k25.lz4.hist',skip=8)
#malbac=read.table('S01.k25.lz4.hist',skip=8)
indata=read.table(infile,skip=8)

getdata <- function(x,range) {
	thefreq <- x[,1]
	cnt <- as.numeric(x[,2])
	theratio <- x[,3]
	thedata <- cbind(thefreq,cnt)
	flag <- thedata[,1] <= range
	#RVAL <- cbind(theratio[flag],thefreq[flag])
	RVAL <- cbind(cnt[flag],thefreq[flag])
	RVAL
}

usedata1 <- getdata(blood,testXrange)
usedata2 <- getdata(indata,testXrange)
#usedata3 <- getdata(malbac,testXrange)

#xx=cbind(usedata2[,1],usedata1[,1])
xx=usedata2[,1];
yy=usedata1[,1]
#yy=cbind(usedata1[,1],usedata2[,1])
#aa=chisq.test(usedata2[,1], p = usedata1[,1], rescale.p = TRUE)
aa=chisq.test(cbind(xx,yy))
#bb=chisq.test(yy)

bb<-wilcox.test(xx,yy,alternative='two',paired=TRUE,conf.int=TRUE)

print(usedata2[,1])
print(usedata1[,1])
#print(xx)
print(c(sum(xx),sum(yy)))
print(aa)
print(bb)
"
scaleres <- myfit(usedata2, par = list(lambda = 3))
scaleres$fitted <- as.vector(usedata1[, 1])
print('===scaledmyfit===')
print(scaleres)
dbg1 <- function(x) {
	arg <- deparse(substitute(x))
	RVAL <- cbind(x$observed,x$fitted,(x$observed-x$fitted)/x$fitted)
	#RVAL <- round(RVAL,digits = 6)
	colnames(RVAL) <- c(paste0(arg,'.ob'),paste0(arg,'.fi'),'(O-E)/E')
	rownames(RVAL) <- 1+ 1:nrow(RVAL)
	RVAL
}
print(cbind(dbg1(scaleres)))
scalesum <- summary(scaleres)
"
mm=pdf2cdf(xx,T,T)
nn=pdf2cdf(yy,T,T)

ksres <- ks.test(xx,yy,alternative='t')
print(ksres)
ksres <- ks.test(mm,nn,alternative='t')
print(ksres)

"
./kmercmp.r blood.fq.k25.gz.hist.fixed S1.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed S2.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed S3.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed S23.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed S24.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed S28.fq.k25.gz.hist.fixed

./kmercmp.r S2.fq.k25.gz.hist.fixed S1.fq.k25.gz.hist.fixed
./kmercmp.r blood.fq.k25.gz.hist.fixed mlbacDonor.k25.lz4.hist.fixed
./kmercmp.r mlbacDonor.k25.lz4.hist.fixed S23.fq.k25.gz.hist.fixed

"
