DATa <- read.table("bamrsplot.tsv",skip=1)
DATb <- read.table("depstat/nrss1m.tsv",skip=3)
for(k in c(1:22)) {
     png(filename = paste0("mixplot.chr",k,".png"),
          width = 3780, height = 2835, units = "px", pointsize = 96)
     par(mar=c(2, 4, 2, 0.5))     # c(bottom, left, top, right)
     layout(rbind(1,2,3), heights=c(7,7,1))
          data <- subset(DATa, V1==k);
          ymax <- 400;
          thelwd = 8
          plot(data[,2],data[,3]/1000,ylim=c(0,ymax),lwd=thelwd,
               main=paste0("Chr",k," Reads distribution"),xlab="M",ylab="Reads Count (k)",type="l",col="red",xaxs="r",yaxs="r")
          lines(data[,2],data[,4]/1000,type="l",col="blue",lwd=thelwd)
          lines(data[,2],data[,5]/1000,type="l",col="green",lwd=thelwd)
          #lines(data[,2],data[,6]/1000,type="l",col="black",lwd=thelwd)
          lines(data[,2],data[,7]/1000,type="l",col="orange",lwd=thelwd)
          lines(data[,2],data[,8]/1000,type="l",col="purple",lwd=thelwd)
          lines(data[,2],data[,9]/1000,type="l",col="gray",lwd=thelwd)
         
          lines(data[,2],data[,6]/1000,type="l",col="black",lwd=thelwd)
          box(lwd=thelwd)
          axis(1,lwd=thelwd)
          axis(2,lwd=thelwd)

          data <- subset(DATb, V1==k);
          X <- data[,24:30]
          YY <- X[cbind(row(X)[which(!X == 0)], col(X)[which(!X == 0)])]
          ValueMax <- mean( YY )
          ymax <- 5*ValueMax
          MinusRatio <- -25*ceiling(ValueMax/10)
          plot(data[,2],(data[,10]*MinusRatio),ylim=c(MinusRatio,ymax), main=paste0("Chr",k," Bases distribution"),xlab="",ylab="UnCovered Rate      Base Coverage Depth",type="l",pch=1 ,col="red",lwd=thelwd,yaxt='n')
          lines(data[,2],data[,24],type="l",col="red",lwd=thelwd)
          lines(data[,2],(data[,11]*MinusRatio),type="l",col="blue",lwd=thelwd)
          lines(data[,2],data[,25],type="l",col="blue",lwd=thelwd)
          lines(data[,2],(data[,12]*MinusRatio),type="l",col="green",lwd=thelwd)
          lines(data[,2],data[,26],type="l",col="green",lwd=thelwd)
          #lines(data[,2],(data[,13]*MinusRatio),type="l",col="black",lwd=thelwd)
          #lines(data[,2],data[,27],type="l",col="black",lwd=thelwd)
          lines(data[,2],(data[,14]*MinusRatio),type="l",col="orange",lwd=thelwd)
          lines(data[,2],data[,28],type="l",col="orange",lwd=thelwd)
          lines(data[,2],(data[,15]*MinusRatio),type="l",col="purple",lwd=thelwd)
          lines(data[,2],data[,29],type="l",col="purple",lwd=thelwd)
          lines(data[,2],(data[,16]*MinusRatio),type="l",col="gray",lwd=thelwd)
          lines(data[,2],data[,30],type="l",col="gray",lwd=thelwd)

          lines(data[,2],(data[,13]*MinusRatio),type="l",col="black",lwd=thelwd)
          lines(data[,2],data[,27],type="l",col="black",lwd=thelwd)

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
     legend('center','groups',c("Sperm23","Sperm24","Sperm28","Donor","SpermS01","SpermS02","SpermS03"), lty=seq(1,1,length=7),lwd=24,col=c("red","blue","green","black","orange","purple","gray"),bty="n",horiz=TRUE)
     dev.off()
}
