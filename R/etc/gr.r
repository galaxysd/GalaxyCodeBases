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
sdat <- new.env()

for(i in types) {
	sdat[[i]] <- t(subset(dat, dat[,3] == i) )
}

# ./gr.r gr.tsv
