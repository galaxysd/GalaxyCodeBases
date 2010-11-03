#!/bin/env Rscript
self='./qtl.R';
argv = commandArgs(T);

if (is.null(argv) | length(argv)<2) {
  cat("Error: No options. Run",self,"<gen_rot.csv> <phe_rot.csv> <output_prefix>.{log,*.png,*.txt,*.dump}\n")
  q(status=1)
}

# argv=c('gen_rot.csv','fv1_phe.csv','test1')

#################################################
outputSO <- function (object, threshold, perms, alpha, lodcolumn = 1, pvalues = FALSE, df = FALSE, ...) {
    if (!any(class(object) == "scanone")) 
        stop("Input should have class \"scanone\".")
    format <- "onepheno"
    ncol.object <- ncol(object) - 2
    cn.object <- colnames(object)[-(1:2)]
    if (format != "onepheno" && !missing(lodcolumn)) 
        warning("lodcolumn ignored except when format=\"onepheno\".")
    if (!missing(perms)) {
        if ("scantwoperm" %in% class(perms)) 
            perms <- scantwoperm2scanoneperm(perms)
        else if (!("scanoneperm" %in% class(perms))) 
            warning("perms need to be in scanoneperm format.")
    }
    if (missing(perms) && !missing(alpha)) 
        stop("If alpha is to be used, permutation results must be provided.")
    if (!missing(threshold) && !missing(alpha)) 
        stop("Only one of threshold and alpha should be specified.")
    if (format == "onepheno") {
        if (!missing(lodcolumn) && length(lodcolumn) > 1) {
            warning("With format=\"onepheno\", lodcolumn should have length 1.")
            lodcolumn <- lodcolumn[1]
        }
        if (lodcolumn < 1 || lodcolumn > ncol.object) 
            stop("lodcolumn should be between 1 and no. LOD columns.")
    }
    if (!missing(alpha) && length(alpha) > 1) {
        warning("alpha should have length 1.")
        alpha <- alpha[1]
    }
    if (!missing(perms)) {
        if ("xchr" %in% names(attributes(perms))) {
            ncol.perms <- ncol(perms$A)
            cn.perms <- colnames(perms$A)
        }
        else {
            ncol.perms <- ncol(perms)
            cn.perms <- colnames(perms)
        }
        if (ncol.object != ncol.perms) {
            if (ncol.perms == 1) {
                origperms <- perms
                if ("xchr" %in% names(attributes(perms))) {
                  for (j in 2:ncol.object) {
                    perms$A <- cbind(perms$A, origperms$A)
                    perms$X <- cbind(perms$X, origperms$X)
                  }
                  cn.perms <- colnames(perms$A) <- colnames(perms$X) <- cn.object
                }
                else {
                  for (j in 2:ncol.object) perms <- cbind(perms, 
                    origperms)
                  cn.perms <- colnames(perms) <- cn.object
                }
                warning("Just one column of permutation results; reusing for all LOD score columns.")
            }
            else {
                if (ncol.object == 1) {
                  warning("Using just the first column in the perms input")
                  if ("xchr" %in% names(attributes(perms))) {
                    perms$A <- perms$A[, 1, drop = FALSE]
                    perms$X <- perms$X[, 1, drop = FALSE]
                  }
                  else {
                    clp <- class(perms)
                    perms <- perms[, 1, drop = FALSE]
                    class(perms) <- clp
                  }
                }
                else stop("scanone input has different number of LOD columns as perms input.")
            }
        }
        if (!all(cn.object == cn.perms)) 
            warning("Column names in scanone input do not match those in perms input.")
    }
    if (missing(perms) && pvalues) {
        warning("Can show p-values only if perms are provided.")
        pvalues <- FALSE
    }
    if (!missing(perms) && "df" %in% names(attributes(object)) && 
        "df" %in% names(attributes(perms))) {
        df.o <- attr(object, "df")
        df.p <- attr(perms, "df")
        if (length(df.o) != length(df.p) || any(df.o != df.p)) 
            warning("Degrees of freedom in input object and perms do not match.")
    }
    chr <- as.character(object[, 1])
    if (format == "onepheno") {
        lodcolumn <- lodcolumn + 2
        wh <- 1:length(object[,1])
        thechr <- as.character(object[wh, 1])
        if (!missing(threshold)) 
            wh <- wh[object[wh, lodcolumn] > threshold]
        else if (!missing(alpha)) {
            thr <- summary(perms, alpha)
            if ("xchr" %in% names(attributes(perms))) {
                thr <- sapply(thr, function(a, b) a[, b], lodcolumn - 
                  2)
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]
                xchr <- thechr %in% xchr
                wh <- wh[(!xchr & object[wh, lodcolumn] > thr[1]) | 
                  (xchr & object[wh, lodcolumn] > thr[2])]
            }
            else {
                thr <- thr[, lodcolumn - 2]
                wh <- wh[object[wh, lodcolumn] > thr]
            }
        }
        result <- object[wh, ]
    }
    if (pvalues && nrow(result) > 0) {
        rn <- rownames(result)
        if ("xchr" %in% names(attributes(perms))) {
            xchr <- attr(perms, "xchr")
            xchr <- names(xchr)[xchr]
            xchr <- as.character(result[, 1]) %in% xchr
            L <- attr(perms, "L")
            Lt <- sum(L)
            pval <- vector("list", ncol.object)
            for (i in 1:ncol.object) {
                if (format == "allpeaks") 
                  thecol <- i * 2 + 1
                else thecol <- i + 2
                pval[[i]] <- rep(0, length(xchr))
                if (any(xchr)) 
                  pval[[i]][xchr] <- sapply(result[xchr, thecol], 
                    function(a, b, rat) 1 - mean(b < a)^rat, 
                    perms$X[, i], Lt/L[2])
                if (any(!xchr)) 
                  pval[[i]][!xchr] <- sapply(result[!xchr, thecol], 
                    function(a, b, rat) 1 - mean(b < a)^rat, 
                    perms$A[, i], Lt/L[1])
            }
        }
        else {
            pval <- vector("list", ncol.object)
            for (i in 1:ncol.object) {
                if (format == "allpeaks") 
                  thecol <- i * 2 + 1
                else thecol <- i + 2
                pval[[i]] <- sapply(result[, thecol], function(a, 
                  b) mean(b >= a), perms[, i])
            }
        }
        if (format == "allpeaks") {
            temp <- as.data.frame(matrix(nrow = nrow(result), 
                ncol = ncol.object * 3 + 1))
            names(temp)[1] <- names(result)[1]
            temp[, 1] <- result[, 1]
            for (i in 1:ncol.object) {
                names(temp)[i * 3 + (-1:1)] <- c(names(result)[i * 
                  2 + (0:1)], "pval")
                temp[, i * 3 - 1:0] <- result[, i * 2 + (0:1)]
                temp[, i * 3 + 1] <- pval[[i]]
            }
        }
        else {
            temp <- as.data.frame(matrix(nrow = nrow(result), 
                ncol = ncol.object * 2 + 2))
            names(temp)[1:2] <- names(result)[1:2]
            temp[, 1:2] <- result[, 1:2]
            for (i in 1:ncol.object) {
                names(temp)[i * 2 + 1:2] <- c(names(result)[i + 
                  2], "pval")
                temp[, i * 2 + 1] <- result[, i + 2]
                temp[, i * 2 + 2] <- pval[[i]]
            }
        }
        result <- temp
        rownames(result) <- rn
    }
    if (df && "df" %in% names(attributes(object))) 
        attr(result, "df") <- attr(object, "df")
    if (!df) 
        attr(result, "df") <- NULL
    if (format == "allpeaks") 
        rownames(result) <- as.character(result$chr)
    class(result) <- c("summary.scanone", "data.frame")
    result
}
#################################################

sink(paste(argv[3],'log',sep='.'));
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
	outperm <- scanone(dat0, n.perm=1000, pheno.col=i, verbose=FALSE)
	print(summary(outperm, alpha=c(0.01, 0.05, 0.20)))
	print(summary(out, perms=outperm, alpha=0.5, pvalues=TRUE))
	a=outputSO(out, perms=outperm, alpha=1, pvalues=TRUE)
	b=outputSO(out, perms=outperm, pvalues=TRUE)
	write.table(a,paste(argv[3],phe[i],'txt',sep='.'))
	write.table(b,paste(argv[3],phe[i],'dump',sep='.'))
	png(paste(argv[3],phe[i],'png',sep='.'),16000,1200)
	plot(out, ylab="LOD score",main='Raw SNP Markers')
	dev.off()

	cat(paste('ScanOne NP for [',phe[i],']:\n',sep=''))
	out.np <- scanone(dat, model="np", pheno.col=i)
	outperm.np <- scanone(dat, model="np", n.perm=1000, pheno.col=i, verbose=FALSE)
	print(summary(outperm.np, alpha=c(0.01, 0.05, 0.20)))
	print(summary(out.np, perms=outperm.np, alpha=0.5, pvalues=TRUE))
	a=outputSO(out.np, perms=outperm.np, alpha=1, pvalues=TRUE)
	b=outputSO(out.np, perms=outperm.np, pvalues=TRUE)
	write.table(a,paste(argv[3],phe[i],'np.txt',sep='.'))
	write.table(b,paste(argv[3],phe[i],'np.dump',sep='.'))
	png(paste(argv[3],phe[i],'np.png',sep='.'),16000,1200)
	plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
	dev.off()
}
cat('Done !\n');
sink();
