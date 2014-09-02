#!/usr/bin/env littler

COLORS <- c("red","orange","#FF1EEA","black","#02CC7B","blue","#3399CC")
SampleIDs <- c("Sperm23","Sperm24","Sperm28","Donor","SpermS01","SpermS02","SpermS03")
COLORS <- c("gray","black","red","orange","#FF1EEA","#02CC7B","blue","#3399CC")
SampleIDs <- c("MDA-Blood","MALBAC-Blood","MDA-Sperm23","MDA-Sperm24","MDA-Sperm28","MALBAC-SpermS01","MALBAC-SpermS02","MALBAC-SpermS03")
#COLORS <- c("red","blue","green","black","orange","purple","gray")
theLen <- length(SampleIDs)

ChrCount <- 22
#ChrCount <- 1

DATa <- read.table("bamrsplot8.tsv.gz",skip=1)
DATb <- read.table("rss1m.tsv.gz",skip=3)
toplot <- function(value=3) {
	for(k in 1:ChrCount) {
#		png(filename = paste0("mixplot.chr",k,".png"),
#		  width = 3780, height = 2835, units = "px", pointsize = 96)
		tiff(filename = paste0("mixplot.chr",k,".tiff"), compression="lzw", 
			width = 2049, height = 1800, units = "px", pointsize = 52)
		par(mar=c(2, 4, 2, 0.5))     # c(bottom, left, top, right)
		layout(rbind(1,2,3), heights=c(7,7,1))
		data <- subset(DATa, V1==k);
		ymax <- 800;
		thelwd = value
		plot(data[,2],data[,3]/1000,ylim=c(0,ymax),lwd=thelwd,
		     main=paste0("Chr",k," Reads distribution"),xlab="M",ylab="Reads Count (k)",type="l",col=COLORS[1],xaxs="r",yaxs="r")

		for (cc in 2:length(SampleIDs) ) {
			lines(data[,2],data[,2+cc]/1000,type="l",col=COLORS[cc],lwd=thelwd)
		}

		box(lwd=thelwd)
		axis(1,lwd=thelwd)
		axis(2,lwd=thelwd)

		chrk <- paste0('chr',k)
		data <- subset(DATb, V1==chrk);
		X <- data[,(3+3*theLen):(2+4*theLen)]
		YY <- X[cbind(row(X)[which(!X == 0)], col(X)[which(!X == 0)])]
		ValueMax <- mean( YY )
		ymax <- 3*ValueMax
		MinusRatio <- -25*ceiling(ValueMax/10)
		plot(data[,2],(data[,3+theLen]*MinusRatio),ylim=c(MinusRatio,ymax), main=paste0("Chr",k," Bases distribution"),xlab="",ylab="UnCovered Rate          Base Coverage Depth",type="l",pch=1 ,col=COLORS[1],lwd=thelwd,yaxt='n')
		lines(data[,2],data[,3+3*theLen],type="l",col=COLORS[1],lwd=thelwd)
		for (cc in 2:length(SampleIDs) ) {
			lines(data[,2],(data[,2+theLen+cc]*MinusRatio),type="l",col=COLORS[cc],lwd=thelwd)
			lines(data[,2],data[,2+3*theLen+cc],type="l",col=COLORS[cc],lwd=thelwd)
		}
		box(lwd=thelwd)
		axis(1,lwd=thelwd)
		#axis(2,lwd=thelwd,labels=NA)
		axis(2,lwd=thelwd,at=axTicks(2)[axTicks(2)>=0])
		#TmpTmp <- c(axTicks(2)[axTicks(2)<0], MinusRatio)
		TmpTmp <- c(.25,.5,.75)*MinusRatio
		axis(2,lwd=thelwd,at=c(TmpTmp,MinusRatio),labels=c(paste0(100*TmpTmp/MinusRatio,'%'),1),las=1 )

		# setup for no margins on the legend
		par(mar=c(0, 0, 0, 0))
		# c(bottom, left, top, right)
		plot.new()
		legend('center','groups',SampleIDs, lty=seq(1,1,length=theLen),lwd=16,col=COLORS,bty="n",horiz=TRUE,cex=0.56)
		dev.off()
	}
}

toplot(3)
