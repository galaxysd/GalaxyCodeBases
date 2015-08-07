#!/usr/bin/env littler

if (!interactive()) {
	if (is.null(argv) | length(argv)<3) {
		cat("Usage: ./vcfdensity.r <WinSize> <input.vcf.gz> <output.prefix>\n")
		q(status=1)
	} else {
		WinSize <- as.integer(argv[1])
		if (is.na(WinSize)) WinSize <- 0L
		InVCFgz <- trimws(argv[2],'both')
		OutP <- trimws(argv[3],'both')
		cat("[!] WinSize=[",WinSize,"], [",InVCFgz,"]->[",OutP,"].*\n",sep='')
		if (!file.exists(InVCFgz)) {
			cat("[x] File not found:[", InVCFgz, "] !\n",sep='')
			q(status=-1)
		}
		if ( WinSize < 1L ) {
			cat("[x] WinSize=[", WinSize, "] is smaller than [1] !\n",sep='')
			q(status=-1)
		}
	}
}

library('data.table')
suppressPackageStartupMessages(library('zoo'))
#tab5rows <- read.table(pipe(paste0("gzip -dc ",InVCFgz," | cut -f1,2")),sep="\t",nrows=5)
#classes <- sapply(tab5rows, class)
classes <- c(V1="factor",V2="integer")
#tabAll <- read.table(pipe("zcat vcf.gz | cut -f1,2|head -300"),sep="\t", colClasses = classes,col.names=c('Chr','Pos'))
tabAll <- fread(paste0("gzip -dc ",InVCFgz,"|awk '!/^#|\tINDEL;/'"),header=F,stringsAsFactors=T,sep="\t",autostart=100,select=c(1,2), colClasses=classes, data.table=T,verbose=F)
# grep -ve '^#' 也可，但 MacOS X 下没有 grep -vP '\tINDEL'
setnames(tabAll,1,'Chr')
setnames(tabAll,2,'Pos')
print(head(tabAll))
cat("...\t...\t...\n")
print(tail(tabAll))

dorolling <- function(x, rollwin) {
	chrdat <- integer(max(x))
	chrdat[x] <- 1L
	#thelen <- length(chrdat)
	#length(chrdat) <- ceiling(thelen/WinSize)*WinSize	# 补齐末端会造成 bias
	#chrdat[is.na(chrdat)] <- 0L
	res0 <- rollapply(chrdat, rollwin, sum, by = rollwin)
	return(res0)
}
resArr <- tabAll[, dorolling(Pos,rollwin=WinSize), by=Chr]
resAll <- unlist(resArr$V1,use.names=F)
cat("[!] Stat done.\n")

tbl <- table(resAll)
print(tbl)
write.table(tbl, paste0(OutP,".tsv"), sep = "\t", quote=F,row.names=F,col.names=F)

reshist <- hist(resAll,plot=F)
pdf(file=paste0(OutP,".pdf"),title='Histogram of VCF Density')
plot(reshist,freq=F,main='Histogram of SNP windowed density',xlab=paste0('SNP Count in every ',WinSize,' bps'))
dev.off()

resPhist <- hist(resAll[resAll!=0],plot=F)
pdf(file=paste0(OutP,".nonZero.pdf"),title='Histogram of VCF Density')
plot(resPhist,freq=F,main='Histogram of SNP windowed density (+)',xlab=paste0('SNP Count in every ',WinSize,' bps'))
dev.off()

