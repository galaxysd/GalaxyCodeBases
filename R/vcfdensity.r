#!/usr/bin/env littler

library('parallel')
NumofCore <- detectCores(logical = TRUE)

# 这是实验室同学委托咱写的，统计给定尺寸的窗口内，SNP出现次数的分布。
# 咱觉得那篇文章统计这个纯属无聊，只能证明该物种基因组很大，能被分成50个以上的100k的窗口。即“中心极限定理”。
# 好吧，刚才确认定理名时看到：“独立同分布的中心极限定理”，于是，算是证明该个体，在100k的尺度上，没有明显的连锁，吧。额，貌似有点用？
# 貌似美元符号会按LaTex解析，额，反正下面有原文本，大家将就看吧。为这个折腾转义也太EP了。
# 话说，SNP频率不独立，和连锁没关系吧。那是证明啥？在这个尺度上没有“鱼”和“饵”的关系？本来真核生物中，也就酵母等微生物才把同一代谢通路的基因挨着放吧。于是，这统计到底有啥用？
# http://pastebin.com/Ks0F3r3y
# https://twitter.com/galaxy001/status/630356226782564353

if (!interactive()) {
	if (is.null(argv) | length(argv)<3) {
		cat("Usage: ./vcfdensity.r <WinSize> <input.vcf.gz> <output.prefix>\n")
		q(status=1)
	} else {
		WinSize <- as.integer(argv[1])
		if (is.na(WinSize)) WinSize <- 0L
		InVCFgz <- trimws(argv[2],'both')
		OutP <- trimws(argv[3],'both')
		cat("[!] WinSize=[",WinSize,"], [",InVCFgz,"]->[",OutP,"].*, Core:[",NumofCore,"]\n",sep='')
		if (!file.exists(InVCFgz)) {
			cat("[x] File not found:[", InVCFgz, "] !\n",sep='')
			q(status=-1)
		}
		if ( WinSize < 1L ) {
			cat("[x] WinSize=[", WinSize, "] is smaller than [1] !\n",sep='')
			q(status=-1)
		}
	}
} else {	# for `source("vcfdensity.r")`
	InVCFgz <- 't.vcf'
	WinSize <- 100L
	OutP <- 'test'
	cat("[!] WinSize=[",WinSize,"], [",InVCFgz,"]->[",OutP,"].*, Core:[",NumofCore,"]\n",sep='')
}

library('data.table')
#suppressPackageStartupMessages(library('zoo'))
# http://www.r-bloggers.com/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
wapply <- function(x, width, by = NULL, FUN = NULL, ...) {
	FUN <- match.fun(FUN)
	if (is.null(by)) by <- width
	lenX <- length(x)
	SEQ1 <- seq(1, lenX - width + 1, by = by)
	SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
	OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
	#OUT <- mclapply(SEQ2, function(a) FUN(x[a], ...), mc.cores = getOption("mc.cores", NumofCore))
	OUT <- base:::simplify2array(OUT, higher = TRUE)
	return(OUT)
}

options(datatable.verbose=T)

#tab5rows <- read.table(pipe(paste0("gzip -dc ",InVCFgz," | cut -f1,2")),sep="\t",nrows=5)
#classes <- sapply(tab5rows, class)
classes <- c(V1="factor",V2="integer")
#tabAll <- read.table(pipe("zcat vcf.gz | cut -f1,2|head -300"),sep="\t", colClasses = classes,col.names=c('Chr','Pos'))

if (grepl("\\.gz$",InVCFgz,ignore.case=T)) {
	fileIn <- paste0("gzip -dc ",InVCFgz,"|awk '!/^#|\tINDEL;/'|cut -f1,2")
# grep -ve '^#' 也可，但 MacOS X 下没有 grep -vP '\tINDEL'
# Piped file will be in `/dev/shm/file46583ba3b517` tempory, thus add `|cut -f1,2` to reduce its footprint.
} else if (grepl("\\.vcf$",InVCFgz,ignore.case=T)) {
	fileIn <- InVCFgz
} else {
	fileIn <- paste0("bcftools view -H -vsnps ",InVCFgz,"|cut -f1,2")
}

tabAll <- fread(fileIn, header=F,stringsAsFactors=T,sep="\t",autostart=100,select=c(1,2), colClasses=classes, data.table=T)
setnames(tabAll,1,'Chr')
setnames(tabAll,2,'Pos')
print(head(tabAll))
cat("...\t...\t...\n")
setkey(tabAll,Chr)
print(tail(tabAll))

dorolling <- function(x, rollwin) {
	chrdat <- integer(0)	# No need to `integer(max(x))` as we use `sum(na.rm=T)`
	chrdat[x] <- 1L
	#thelen <- length(chrdat)
	#length(chrdat) <- ceiling(thelen/WinSize)*WinSize	# 补齐末端会造成 bias
	#chrdat[is.na(chrdat)] <- 0L
	res0 <- wapply(chrdat, rollwin, FUN = sum,na.rm = TRUE)
	return(res0)
}
#resArr <- tabAll[, dorolling(Pos,rollwin=WinSize), by=Chr]
#resAll <- unlist(resArr$V1,use.names=F)

#resArr <- tabAll[, .(Density=dorolling(Pos,rollwin=WinSize)), by=Chr]
# lapply optimization changed j from 'lapply(.SD, dorolling, rollwin = WinSize)' to 'list(dorolling(Pos, rollwin = WinSize))'

# http://stackoverflow.com/questions/19082794/speed-up-data-table-group-by-using-multiple-cores-and-parallel-programming
resArr <- mclapply(unique(tabAll[["Chr"]]),
	function(x) {
		tabAll[.(x),.(Density=dorolling(Pos,rollwin=WinSize)), by=Chr]
	}, mc.cores=NumofCore)
resArr <- rbindlist(resArr)

resAll <- unlist(resArr$Density,use.names=F)
cat("[!] Stat done.\n")

tbl <- table(resAll)
print(head(tbl,80))
write.table(tbl, paste0(OutP,".tsv"), sep = "\t", quote=F,row.names=F,col.names=F)

reshist <- hist(resAll,plot=F)
pdf(file=paste0(OutP,".pdf"),title='Histogram of VCF Density')
plot(reshist,freq=F,main='Histogram of SNP windowed density',xlab=paste0('SNP Count in every ',WinSize,' bps'))
dev.off()

resPhist <- hist(resAll[resAll!=0],plot=F)
pdf(file=paste0(OutP,".nonZero.pdf"),title='Histogram of VCF Density')
plot(resPhist,freq=F,main='Histogram of SNP windowed density (+)',xlab=paste0('SNP Count in every ',WinSize,' bps'))
dev.off()

