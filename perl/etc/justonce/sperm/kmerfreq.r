#!/usr/bin/env r
if (is.null(argv) | length(argv)<2) {
  cat("Usage: kmerfreq.r hist_file coverage\n")
  q()
}
testXrange <- as.numeric(argv[1])
infile <- argv[2]
avgcvg <- as.numeric(argv[3])

if(!exists("D__myfit"))
	source("myfit.R")

#inid <- unlist(strsplit(infile,'.',TRUE))[1]
inid <- strsplit(infile,'.',TRUE)[[1]][1]

maxXrange <- 25
"
testXrange <- 15
infile <- 'blood.fq.k25.gz.hist.fixed'
avgcvg <- 3.27885
"
library('vcd')

aa=read.table(infile,skip=8)
thefreq <- aa[,1]
cnt <- as.numeric(aa[,2])
theratio <- aa[,3]

thedata <- cbind(thefreq,cnt)
flag <- thedata[,1] <= testXrange
#usedata <- cbind(cnt[flag],thefreq[flag])
usedata <- cbind(theratio[flag],thefreq[flag])

fitres <- goodfit(usedata,type='poisson', par = list(lambda = avgcvg))
#fitres <- goodfit(usedata,type='poisson', method = 'MinChisq')

#sumres <- summary(fitres)

print.myfit <- function(x, ...)
{
    cat(paste("\nObserved and fitted values for", x$type, "distribution\n"))
    if(x$method[1] == "fixed")
      cat("with fixed parameters \n\n")
    else
      cat(paste("with parameters estimated by `", paste(x$method, collapse = " "), "' \n\n", sep = ""))
    RVAL <- cbind(x$count, x$observed, x$fitted)
    colnames(RVAL) <- c("count", "observed", "fitted")
    rownames(RVAL) <- rep("", nrow(RVAL))
    print(RVAL)
    invisible(x)
}
fitted.myfit <- function(object, ...)
{
  object$fitted
}
summary.myfit <- function(object, ...)
{
    df <- object$df
    obsrvd <- object$observed
    count <- object$count
    expctd <- fitted(object)

    G2 <- sum(ifelse(obsrvd == 0, 0, obsrvd * log(obsrvd/expctd))) * 2

    n <- length(obsrvd)
    switch(object$type,
    "poisson" = { pfun <- "ppois" },
    "binomial" = { pfun <- "pbinom" },
    "nbinomial" = { pfun <- "pnbinom" })
    p.hat <- diff(c(0, do.call(pfun, c(list(q = count[-n]), object$par)), 1))
    #expctd <- p.hat * sum(obsrvd)
    X2 <- sum((obsrvd - expctd)^2/expctd)

    names(G2) <- "Likelihood Ratio"
    names(X2) <- "Pearson"
    if(any(expctd < 5) & object$method[1] != "ML")
        warning("Chi-squared approximation may be incorrect")

    RVAL <- switch(object$method[1],
      "ML" = { G2 },
      "MinChisq" = { X2 },
      "fixed" = { c(X2, G2) }
    )

    RVAL <- cbind(RVAL, df, pchisq(RVAL, df = df, lower.tail = FALSE))
    colnames(RVAL) <- c("X^2", "df", "P(> X^2)")

    cat(paste("\n\t Goodness-of-fit test for", object$type, "distribution\n\n"))
    print(RVAL)
    invisible(RVAL)
}

myfit <- function(x, method = c("ML", "MinChisq"), par = NULL) {
	if (is.vector(x)) {
		x <- table(x)
	}
	if (is.table(x)) {
		if (length(dim(x)) > 1)
			stop("x must be a 1-way table")
		freq <- as.vector(x)
		count <- as.numeric(names(x))
	} else {
		if (!(!is.null(ncol(x)) && ncol(x) == 2))
			stop("x must be a 2-column matrix or data.frame")
		freq <- as.vector(x[, 1])
		count <- as.vector(x[, 2])
	}
	#ncount <- 0:max(count)
	n <- length(count)
	df <- -1
	method <- match.arg(method)
	if (!is.null(par)) {
		if (!is.list(par)) stop("`par' must be a named list")
		if (!(names(par) == "lambda")) stop("`par' must specify `lambda'")
		par <- par$lambda
		method <- "fixed"
	} else if (method == "ML") {
		df <- df - 1
		par <- weighted.mean(count, freq)
	} else if (method == "MinChisq") {
		df <- df - 1
		chi2 <- function(x) {
			p.hat <- diff(c(0, ppois(count[-n], lambda = x),
			  1))
			expected <- sum(freq) * p.hat
			sum((freq - expected)^2/expected)
		}
		par <- optimize(chi2, range(count))$minimum
	}
	par <- list(lambda = par)
	p.hat <- dpois(count, lambda = par$lambda)
	#expected <- sum(freq)/sum(p.hat) * p.hat
	#expected <- p.hat
	zoomratio <- sum(freq) / sum(p.hat)
	expected <- p.hat * zoomratio
	print(c(sum(freq),sum(p.hat),zoomratio))
	print(cbind(p.hat,expected))
	#relsum <- sum(expected)/sum(freq);
	#print(relsum);
	#freq <- freq*relsum;
	df <- switch(method[1], MinChisq = {
		length(freq) + df
	}, ML = {
		sum(freq > 0) + df
	}, fixed = {
		c(length(freq), sum(freq > 0)) + df
	})
	RVAL <- list(observed = freq, count = count, fitted = expected,
		type = "poisson", method = method, df = df, par = par, zoom=zoomratio)
	class(RVAL) <- "myfit"
	RVAL
}

selfres <- myfit(usedata, par = list(lambda = avgcvg))
print('===myfit===')
print(selfres)
selfsum <- summary(selfres)
#print(selfsum)
selfchi <- chisq.test(cbind(selfres$observed,selfres$fitted))
print(selfchi)

newres<-list(observed = fitres$observed[c(-1,-2)], count = fitres$count[c(-1,-2)], fitted = fitres$fitted[c(-1,-2)],
		type = "poisson", method = fitres$method, df = fitres$df - 2, par = fitres$par)
class(newres) <- "goodfit"
print('===goodfit===')
print(newres)
sumres <- summary(newres)
#print(sumres)

print(cbind(newres$observed,newres$fitted,(newres$observed-newres$fitted)/newres$fitted, selfres$fitted,(selfres$observed-selfres$fitted)/selfres$fitted))

tiff(filename = paste0(infile,".tiff"), compression="lzw", width = 683, height = 683)#, units = "px", pointsize = 52)
par(mar=c(5, 6, 2, 2),ps=20,family='sans')

yy=aa[,3]
xx=aa[,1]
#xx <- usedata[,2]	# count
#yy <- usedata[,1]	# freq

mean2=avgcvg
xxx=seq(0,max(xx))
zz=dpois(xxx,mean2);
maxY=max(zz,yy,zz*selfres$zoom)
maxY<-0.32
#lines(xx,zz,xlim=c(0,60),type='l');

#barplot(yy,xlim=c(0,20),ylim=c(0,0.212),xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
plot(xx,yy,type='h',lwd=20,xlim=c(1,maxXrange),ylim=c(0,maxY),
	main = paste0("Pearson goodness of fit, p=",selfsum[1,3]),
	xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
# dpois(3,3.62)=0.2117524
legend("topright",pch=-1,lty=1,col=c("navy",'green3',"red"), x.intersp = 1, y.intersp = 2,cex=1,lwd=c(11,7,5),
	legend= c(inid,paste0("lamda=",avgcvg),paste0("scaledby=",selfres$zoom)) )
axis(1, at = seq(0, maxXrange, by = 5),lwd=3,cex=1)
#lines(fitres$count,fitres$observed,xlim=c(0,60),type='l',col='blue')
#lines(fitres$count,fitres$fitted,xlim=c(0,60),type='l',col='red',lwd=6)
lines(xxx,zz,xlim=c(0,maxXrange),type='l',col=rgb(0,0.804,0,0.8),lwd=7)
lines(selfres$count,selfres$fitted,xlim=c(0,60),type='l',col=rgb(1,0,0,0.75),lwd=5)
lines(xxx,zz*selfres$zoom,xlim=c(0,maxXrange),type='l',col=rgb(1,0,0,0.8),lwd=3)
dev.off()

print(selfsum)

"
#./kmerfreq.r mlbacDonor.k25.lz4.hist 6.85
#./kmerfreq.r mdaS23.k25.lz4.hist 4.34
#./kmerfreq.r S01.k25.lz4.hist 3.62

#                         X^2 df  P(> X^2)
#Pearson           1.56832025 18 0.9999998
#Likelihood Ratio -0.06964308 18 1.0000000

#                          X^2 df P(> X^2)
#Pearson          7285.4140548 18        0
#Likelihood Ratio    0.3987449 18        1

#                          X^2 df  P(> X^2)
#Pearson          2.511278e+05 18 0.0000000
#Likelihood Ratio 2.540676e+00 18 0.9999924

./kmerfreq.r 15 S1.fq.k25.gz.hist.fixed 3.29571
./kmerfreq.r 15 S2.fq.k25.gz.hist.fixed 3.28174
./kmerfreq.r 15 S3.fq.k25.gz.hist.fixed 2.94459
./kmerfreq.r 15 S23.fq.k25.gz.hist.fixed 3.26044
./kmerfreq.r 15 S24.fq.k25.gz.hist.fixed 3.27885
./kmerfreq.r 15 S28.fq.k25.gz.hist.fixed 3.30884
./kmerfreq.r 15 blood.fq.k25.gz.hist.fixed 3.27885
"
