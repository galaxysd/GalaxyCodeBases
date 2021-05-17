#!/usr/bin/env littler

require(readxl)

# 设置数值显示位数
options(scipen = 200)

argv <- if (!exists("argv")) commandArgs(TRUE) else as.character(argv)
fpath <- argv[1]
fpattern <- "(NFO|HLA|SNP|MED).xlsx?$"

readXLS <- function(xlsname) {
	read <- read_excel(xlsname,col_types="text",col_names=c('k','v'))
	print(read)
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

