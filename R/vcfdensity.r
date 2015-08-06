#!/usr/bin/env littler

if (!interactive()) {
	if (exists("argv")) {
		if (!is.null(argv) && length(argv)>0) {
			print(argv)
		}
	}
}

library('data.table')
library('zoo')
tab5rows <- read.table(pipe("gzip -dc vcf.gz | cut -f1,2"),sep="\t",nrows=5)
classes <- sapply(tab5rows, class)
#tabAll <- read.table(pipe("zcat vcf.gz | cut -f1,2|head -300"),sep="\t", colClasses = classes,col.names=c('Chr','Pos'))
tabAll <- fread("gzip -dc vcf.gz|grep -ve '^#' |head -n500000000",header=F,stringsAsFactors=T,sep="\t",autostart=100,select=c(1,2), colClasses=classes, data.table=T,verbose=T)
#tabAll <- fread("cat t.vcf|grep -ve '^#' |head -n500000000",header=F,verbose=T,sep="\t",autostart=100,select=c(1,2),stringsAsFactors=T, colClasses=classes,data.table=T)

print(tabAll)

Poses <- split(tabAll$V2,tabAll$V1)
print(Poses[2])

WinSize <- 1000

resAll <- integer(0)
for (p in Poses) {
	chrdat <- integer(0)
	for (i in p) chrdat[i] <- 1L
	thelen <- length(chrdat)
	length(chrdat) <- ceiling(thelen/WinSize)*WinSize
	chrdat[is.na(chrdat)] <- 0L
	res0 <- rollapply(chrdat, WinSize, sum, by = WinSize)
	resAll <- append(resAll,res0)
}
res <- hist(resAll,plot=F)
tbl <- table(resAll)

print(resAll)
print(res)
print(tbl)





#t=sapply(Poses,mean)

