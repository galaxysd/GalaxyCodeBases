#!/usr/bin/env littler

if (is.null(argv) | length(argv)<1) {
	cat("Usage: ./gr.r <datafile>\n")
	q()
}

#print(argv)
#print(argv[1])
#dat = read.delim('gr.tsv',header=F)

dat = read.delim(argv[1],header=F)
types = levels(dat[,3])
colnames(dat) <- c('Gene','ChrID','Type','Range',seq(1,51))
datacol = colnames(dat)[-1:-4]
sdat <- new.env()
pdf('gr.pdf')

for(i in types) {
	sdat[[i]] <- subset(dat[datacol], dat[,3] == i)
  boxplot(sdat[[i]],use.cols=T,outline=F,varwidth=T,xlab='Zones',ylab='Meth %',main=i)
}

dev.off()
# ./gr.r gr.tsv
