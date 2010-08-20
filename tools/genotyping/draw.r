#!/opt/blc/genome/biosoft/R/bin/Rscript

self='./draw.r';
argv = commandArgs(T);

if (is.null(argv) | length(argv)<1) {
  cat("Error: No options. Run",self,"<Prefix>.{<bin,<border >png} [output_prefix].png\n")
  q(status=1)
}
if (length(argv)==1) outPrefix=argv[1]

geno.bin=read.table(paste(argv[1],'.bin',sep=''))

theborder=as.matrix(read.table(paste(argv[1],'.border',sep='')))

geno.colors <- as.matrix(geno.bin);
geno.colors[is.na(geno.colors)] <- rgb(0,0,0)
geno.colors[geno.colors==0] <- rgb(0,0.7,0)
geno.colors[geno.colors==1] <- rgb(0,0,0.7)
geno.colors[geno.colors==0.5] <- rgb(0.6,0,0)
rils=matrix(theborder,nrow=nrow(geno.bin),ncol=ncol(geno.bin))	# 1:nrow(geno.bin) or theBin$border
poses=matrix(rep(1:ncol(geno.bin),each=nrow(geno.bin)),nrow=nrow(geno.bin),ncol=ncol(geno.bin))
png(paste(outPrefix,'.png',sep=''), 1200, 4800, bg='white')
plot(poses, rils, col=geno.colors,pch=20, ylab='',xlab="RIL index",mar=c(1.1,2.1,2.1,1.1))

for(j in 1:ncol(geno.colors)) {
	for(i in 1:(nrow(geno.colors)-1)) {
		if(geno.colors[i,j]==geno.colors[i+1,j])
			segments(j,rils[i,j],j,rils[i+1,j],col=geno.colors[i,j],lwd=2)
	}
}

for(i in 1:nrow(rils)) {
	segments(1,rils[i,1],ncol(rils),rils[i,ncol(rils)],lwd=1,col=rgb(0,0,0.2,0.4))
}

dev.off()

warnings()

