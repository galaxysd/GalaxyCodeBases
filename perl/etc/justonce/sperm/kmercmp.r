#!/usr/bin/env r
if (is.null(argv) | length(argv)<1) {
  cat("Usage: kmerfreq.r hist_file\n")
  q()
}
infile <- argv[1]

maxXrange <- 30
testXrange <- 20
#infile <- 'mlbacDonor.k25.lz4.hist'
#avgcvg <- 6.85

#library('vcd')

blood=read.table('blood.fq.k25.gz.hist',skip=8)
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

xx=cbind(usedata2[,1],usedata1[,1])
#yy=cbind(usedata1[,1],usedata2[,1])
#aa=chisq.test(usedata2[,1], p = usedata1[,1], rescale.p = TRUE)
aa=chisq.test(xx)
#bb=chisq.test(yy)

print(usedata2[,1])
print(usedata1[,1])
print(xx)
print(c(sum(xx[,1]),sum(xx[,2])))
print(aa)
#print(bb)

#./kmercmp.r mdaS23.k25.lz4.hist
#./kmercmp.r S01.k25.lz4.hist
