#!/bin/env Rscript
self='./qtl.R';
argv = commandArgs(T);

if (is.null(argv) | length(argv)<2) {
  cat("Error: No options. Run",self,"<gen_rot.csv> <phe_rot.csv> <output_prefix>.{logp,*.png}\n")
  q(status=1)
}

# argv=c('gen_rot.csv','fv1_phe.csv','test1')

sink(paste(argv[3],'logp',sep='.'));
library('qtl');
dat=read.cross('csvsr','.',argv[1],argv[2]);

phe=colnames(dat$pheno);
phe=phe[-length(phe)];

for(i in 1:length(phe)) {
  png(paste(argv[3],phe[i],'in.png',sep='.'),640,480);
  plot.pheno(dat,i,col='gray');
  dev.off();
}

png(paste(argv[3],'in.pairs.png',sep='.'),1000,1000);
pairs(jitter( as.matrix(dat$pheno[,1:length(phe)]) ), cex=0.6, las=1)
dev.off();

dat0 <- dat
dat0 <- calc.genoprob(dat0, step=1, error.prob=0.001)
dat <- calc.genoprob(dat, step=1, error.prob=0.01)

for(i in 1:length(phe)) {
	cat(paste('ScanOne Default for [',phe[i],']:\n',sep=''))
	out <- scanone(dat0, pheno.col=i)
	png(paste(argv[3],phe[i],'png',sep='.'),16000,1200)
	plot(out, ylab="LOD score",main='Raw SNP Markers')
	dev.off()

	cat(paste('ScanOne NP for [',phe[i],']:\n',sep=''))
	out.np <- scanone(dat, model="np", pheno.col=i)
	png(paste(argv[3],phe[i],'np.png',sep='.'),16000,1200)
	plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
	dev.off()
}
cat('Done !\n');
sink();
