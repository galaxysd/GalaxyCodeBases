#!/usr/bin/env littler

# R --save < lm.R --args bamrsplot.tsv tsv/nrss1m.tsv read.lm uncover.lm base.lm

# gzcat bamrsplot8.tsv.gz|perl -lane 'print $_ unless $F[0] =~ /[XYM]/i' > bamrsplot8.auto.tsv
# gzcat rss1m.tsv.gz|perl -lane 'print $_ unless $F[0] =~ /[XYM]/i' > rss1m.auto.tsv

names <- c("MDA-Blood","MALBAC-Blood","MDA-Sperm23","MDA-Sperm24","MDA-Sperm28","MALBAC-SpermS01","MALBAC-SpermS02","MALBAC-SpermS03")
theLen <- length(names)

DATa <- read.table("bamrsplot8.auto.tsv",skip=1)
DATb <- read.table("rss1m.auto.tsv",skip=3)

for (i in seq(3,1+theLen)) {
	for (j in seq(i+1,2+theLen)) {
		ta <- cor.test(DATa[,i],DATa[,j],method='pearson')
		tb <- cor.test(DATa[,i],DATa[,j],method='kendall')
		
		cat(as.character(names[i-2]), as.character(names[j-2]), as.character(ta$p.value), as.character(ta$estimate), as.character(tb$p.value), as.character(tb$estimate),"\n",file='read.lm',append=T,sep="\t")
	}
}

for (i in seq(3,2+theLen)){
	var=sd(DATa[,i])/mean(DATa[,i])
	cat(as.character(names[i-2]), as.character(var),"\n",file='read.lm',append=T,sep="\t")		
}

for (i in seq(3+theLen,1+2*theLen)) {
	for (j in seq(i+1,2+2*theLen)) {
		ta <- cor.test(DATb[,i],DATb[,j],method='pearson')
		tb <- cor.test(DATb[,i],DATb[,j],method='kendall')
		cat(as.character(names[i-2-theLen]), as.character(names[j-2-theLen]), as.character(ta$p.value), as.character(ta$estimate), as.character(tb$p.value), as.character(tb$estimate),"\n",file='uncover.lm',append=T,sep="\t")	
	}
	
}
for(i in seq(3+theLen,2+2*theLen)){
	var=sd(DATb[,i])/mean(DATb[,i])
        cat(as.character(names[i-2-theLen]), as.character(var),"\n",file='uncover.lm',append=T,sep="\t")
}

for (i in seq(3+3*theLen,1+4*theLen)) {
	for (j in seq(i+1,2+4*theLen)) {
		ta <- cor.test(DATb[,i],DATb[,j],method='pearson')
		tb <- cor.test(DATb[,i],DATb[,j],method='kendall')
		cat(as.character(names[i-2-3*theLen]), as.character(names[j-2-3*theLen]), as.character(ta$p.value), as.character(ta$estimate), as.character(tb$p.value), as.character(tb$estimate),"\n",file='base.lm',append=T,sep="\t")
	}
}

for(i in seq(3+3*theLen,2+4*theLen)){
        var=sd(DATb[,i])/mean(DATb[,i])
        cat(as.character(names[i-2-3*theLen]),as.character(var),"\n",file='base.lm',append=T,sep="\t")
}
