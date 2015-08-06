#!/usr/bin/env littler

library('data.table')
library('zoo')
tab5rows <- read.table(pipe("zcat vcf.gz | cut -f1,2"),sep="\t",nrows=5)
classes <- sapply(tab5rows, class)
#tabAll <- read.table(pipe("zcat vcf.gz | cut -f1,2|head -300"),sep="\t", colClasses = classes,col.names=c('Chr','Pos'))
#tabAll <- fread("zcat vcf.gz|grep -vP '^#' |head -n500000000",header=F,verbose=T,sep="\t",autostart=100,select=c(1,2),stringsAsFactors=T, colClasses=classes,data.table=T)
tabAll <- fread("cat t.vcf|grep -vP '^#' |head -n500000000",header=F,verbose=T,sep="\t",autostart=100,select=c(1,2),stringsAsFactors=T, colClasses=classes,data.table=T)

print(tabAll)

Poses <- split(tabAll$V2,tabAll$V1)
print(Poses[2])

t=sapply(Poses,mean)
print(t)

WinSize <- 10

a <- c()
a[1]=1
a[10]=1
a[20]=1
a[56]=1
a[58]=1
a[59]=1
a[67]=1
a[is.na(a)] <- 0
t=rollapply(a, WinSize, sum, by = WinSize)
print(a)
print(t)
