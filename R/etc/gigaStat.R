library(RCurl)
library(XML)
library(hash)
library(RColorBrewer)

hyear <- hash()
htype <- hash()
D <- read.table("giga.tsv",sep="\t",quote="",header=TRUE)
nrec <- length(D[,c("DOI")])
vrec <- 0
for(i in 1:nrec) {
	print(sprintf("Processing %d out of %d.",i,nrec))
	doi <- D[i,c("DOI")]
	if(doi == "") {
		next
	}
	site <- "https://academic.oup.com/gigascience/article-lookup/doi/"
	loca <- paste(site, D[i,c("DOI")], sep="")
	header <- c("User-Agent"="Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9.1.6) ","Accept"="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8","Accept-Language"="en-us","Connection"="keep-alive","Accept-Charset"="GB2312,utf-8;q=0.7,*;q=0.7")
	temp <- getURL(loca,httpheader=header,followlocation=T,encoding="UTF-8",ssl.verifypeer = FALSE)
	pagetree <- htmlTreeParse(temp,encoding="UTF-8", useInternalNodes = TRUE,trim=TRUE)
	node <- getNodeSet(pagetree, "//div[@class='citation-date']/text()")
	if(is.null(node)) {
		next
	}
	vrec <- vrec + 1
	ptime <- sapply(node,xmlValue) 
	stime <- strsplit(ptime," ")
	ptime <- stime[[1]][3]
	node <- getNodeSet(pagetree, "//div[@class='article-metadata-tocSections']/a/text()")
	ptype <- sapply(node,xmlValue) 
	node <- getNodeSet(pagetree, path="//meta[@name='citation_author_institution']")
	pinst <- sapply(node,xmlGetAttr,"content")
	ninst <- length(pinst)
	linst <- 0
	for(j in 1:ninst) {
		if(length(grep('BGI', pinst[j])) > 0) {
			linst <- 1
			next
		}
	}

	if(is.null(hyear[[ptime]])) {
		.set(hyear,keys=ptime,values=1)
	} else {
		.set(hyear,keys=ptime,values=hyear[[ptime]]+1)
	}
	
	if(is.null(htype[[ptype]])) {
		.set(htype,keys=ptype,values=1)
	} else {
		.set(htype,keys=ptype,values=htype[[ptype]]+1)
	}
	
	dfdoi <- data.frame(doi)
	dftime <- data.frame(ptime)
	dftype <- data.frame(ptype)
	dfinst <- data.frame(linst)
	dflocal <- merge(dfdoi, dftime)
	dflocal <- merge(dflocal, dftype)
	dflocal <- merge(dflocal, dfinst)
	
	if(i == 1) {
		dfglobal <- dflocal
	} else {
		dfglobal <- rbind(dfglobal, dflocal)
	}
	
	Sys.sleep(1)
}

pdf("Rst.pdf")
kyear <- keys(hyear)
nyear <- as.numeric(kyear)
nyear <- sort(nyear)
vdim <- c(1)
for(i in 1:length(nyear)) {
	vdim[i] <- hyear[[as.character(nyear[i])]]
}
names(vdim) <- as.character(nyear)
barplot(vdim,xlab="Year",ylab="Number",col=brewer.pal(length(vdim), "Paired"),main="Paper number per year")

htime <- hash()
for(i in 1:length(nyear)) {
	.set(htime,keys=as.character(nyear[i]),values=i)
}

par(xpd=NA,mar=par()$mar+c(0, 0, 0, 10))
ktype <- keys(htype)
ktype <- sort(ktype)
mtype <- matrix(0, length(ktype), length(nyear))
rows <- ktype
cols <- as.character(nyear)
dimnames(mtype) <- list(rows,cols)
for(i in 1:length(ktype)) {
	for(j in 1:vrec) {
		if(dfglobal[j,'ptype'] == ktype[i]) {
			tmp <- as.character(dfglobal[j,'ptime'])
			mtype[i,htime[[tmp]]] <- mtype[i,htime[[tmp]]] + 1
		}
	}
}

if(nrow(mtype) >= 3) {
	color <- brewer.pal(nrow(mtype),"Paired")
} else if(nrow(mtype) == 2) {
	color <- c("#A6CEE3","#1F78B4")
} else {
	color <- c("#A6CEE3")
}
barplot(mtype,xlab="Year",ylab="Number",col=color,main="Paper for type per type",beside=TRUE)
legend(par()$usr[2],mean(par()$usr[3:4]),legend=rownames(mtype),fill=color,xpd=NA,xjust=0,yjust=0.5)

minst <- matrix(0, 2, length(nyear))
rows <- c("non-BGI","BGI")
cols <- as.character(nyear)
dimnames(minst) <- list(rows,cols)
for(i in 1:vrec) {
	if(dfglobal[i,'linst'] == 0) {
		tmp <- as.character(dfglobal[i,'ptime'])
		minst[1,htime[[tmp]]] <- minst[1,htime[[tmp]]] + 1
	}
}
for(i in 1:vrec) {
	if(dfglobal[i,'linst'] == 1) {
		tmp <- as.character(dfglobal[i,'ptime'])
		minst[2,htime[[tmp]]] <- minst[2,htime[[tmp]]] + 1
	}
}
if(nrow(minst) >= 3) {
	color <- brewer.pal(nrow(mtype),"Paired")
} else if(nrow(minst) == 2) {
	color <- c("#A6CEE3","#1F78B4")
} else {
	color <- c("#A6CEE3")
}
barplot(minst,xlab="Year",ylab="Number",col=color,main="Paper for BGI/non-BGI per type",beside=TRUE)
legend(par()$usr[2],mean(par()$usr[3:4]),legend=rownames(minst),fill=color,xpd=NA,xjust=0,yjust=0.5)

dev.off()


