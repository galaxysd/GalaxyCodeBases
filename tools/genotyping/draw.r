#!/opt/blc/genome/biosoft/R/bin/Rscript

self='./draw.r';
argv = commandArgs(T);

if (is.null(argv) | length(argv)<1) {
  cat("Error: No options. Run",self,"<Prefix>.{<bin,<rgt >png} [output_prefix].png\n")
  q(status=1)
}
if (length(argv)==1) outPrefix=argv[1] else outPrefix=argv[2]

geno.bin=read.table(paste(argv[1],'.bin',sep=''))
geno.rgt=read.table(paste(argv[1],'.rgt',sep=''))

#theborder=as.matrix(read.table(paste(argv[1],'.border',sep='')))
theborder=as.numeric(rownames(geno.bin))
theborderR=as.numeric(rownames(geno.rgt))

geno.colors <- as.matrix(geno.bin);
geno.colorsR <- as.matrix(geno.rgt);

# A:(0,179,0) B:(0,51,179) AB:(230,230,0) NA:(205,0,0)
# A:(0,0.7,0) B:(0,0.2,0.7) AB:(0.9,0.9,0) NA:(0.8,0,0)

geno.colors[is.na(geno.colors)] <- rgb(0.8,0,0)
geno.colors[geno.colors==0]=rgb(0,0.7,0)
geno.colors[geno.colors==1]=rgb(0,0.2,0.7)
geno.colors[geno.colors==0.5]=rgb(0.9,0.9,0)

Alpha=0.5
AlphaM=255*Alpha
geno.colorsR[is.na(geno.colorsR)] <- rgb(0.8,0,0,Alpha)
geno.colorsR[geno.colorsR==0]=rgb(0,0.7,0,Alpha)
geno.colorsR[geno.colorsR==1]=rgb(0,0.2,0.7,Alpha)
geno.colorsR[geno.colorsR==0.5]=rgb(0.9,0.9,0,Alpha)

rils=matrix(theborder,nrow=nrow(geno.bin),ncol=ncol(geno.bin))	# 1:nrow(geno.bin) or theBin$border
poses=matrix(rep(1:ncol(geno.bin),each=nrow(geno.bin)),nrow=nrow(geno.bin),ncol=ncol(geno.bin))

rilsR=matrix(theborderR,nrow=nrow(geno.rgt),ncol=ncol(geno.rgt))
posesR=matrix(rep(1:ncol(geno.rgt),each=nrow(geno.rgt)),nrow=nrow(geno.rgt),ncol=ncol(geno.rgt))

theX=min(1200,180+ncol(geno.rgt)*16)
theY=min(3600,180+nrow(geno.rgt)*32)

png(paste(outPrefix,'M.png',sep=''), theX, theY, bg='white')
plot(poses, rils, col=geno.colors, pch=19, ps=9, ylab='',xlab="RIL index",mar=c(1.1,2.1,2.1,1.1),ylim=c(1,47283185))
points(posesR, rilsR, col=geno.colorsR,pch='*',ps=5)

for(j in 1:ncol(geno.colors)) {
	for(i in 1:(nrow(geno.colors)-1)) {
		if(geno.colors[i,j]==geno.colors[i+1,j])
			segments(j,rils[i,j],j,rils[i+1,j],col=rgb(t(col2rgb(geno.colors[i,j])),alpha=AlphaM,maxColorValue=255),lwd=2)
	}
}
#rgb(t(col2rgb(geno.colors[i,j])),alpha=AlphaM,maxColorValue=255)
#geno.colors[i,j]

dev.off()
warnings()

