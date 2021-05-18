#!/usr/bin/env littler

require(readxl)
require(tibble)

# 设置数值显示位数
options(scipen = 200)

argv <- if (!exists("argv")) commandArgs(TRUE) else as.character(argv)
fpath <- argv[1]
fpattern <- "(NFO|HLA|SNP|MED).xlsx?$"

readXLS <- function(xlsname) {
	read <- read_excel(xlsname,col_types="text",col_names=FALSE,.name_repair = ~ LETTERS[seq_along(.x)])
	tres <- regexpr('(?<Type>\\w{3}).xls',xlsname,perl=TRUE)
	st <- attr(tres, "capture.start")[1, ]
	tstr <- toupper(substring(xlsname,st,st + attr(tres, "capture.length")[1, ] -1))
	rdat <- as.data.frame(t(column_to_rownames(read,var='A')),row.names=make.names(rep_len(tstr,length(read)-1),unique=TRUE))
	return(rdat)
}

analyze <- function(dirname) {
	#print(c('Fa',dirname))
	fnames <- list.files(file.path(fpath,dirname),pattern=fpattern,ignore.case=TRUE,full.names=TRUE)
	NFOname <- fnames[grep('NFO.xls',fnames,ignore.case=TRUE)[1]]
	HLAname <- fnames[grep('HLA.xls',fnames,ignore.case=TRUE)[1]]
	SNPname <- fnames[grep('SNP.xls',fnames,ignore.case=TRUE)[1]]
	MEDname <- fnames[grep('MED.xls',fnames,ignore.case=TRUE)[1]]
	fnGot <- paste0(dirname,':N=',NFOname,',S=',SNPname,',H=',HLAname,',M=',MEDname)
	readXLS(NFOname)
	readXLS(HLAname)
	readXLS(SNPname)
	reads <- sapply(c(NFOname,HLAname,SNPname),readXLS,USE.NAMES=F)
	print(reads)
	return(c(fnGot,'---'))
}

if (dir.exists(fpath)) {
	flist <- list.dirs(fpath,full.names=FALSE,recursive=FALSE)
	#print(flist)
	ret <- sapply(flist,analyze)
	print(ret)
} else {
	cat(paste("[!]Usage: ./do.r <path>\n"))
}

