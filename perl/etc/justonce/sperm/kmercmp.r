#!/usr/bin/env r
if (is.null(argv) | length(argv)<1) {
  cat("Usage: kmerfreq.r hist_file\n")
  q()
}
infile <- argv[1]

maxXrange <- 30
testXrange <- 20
#infile <- 'mlbacDonor.k25.lz4.hist'
avgcvg <- 6.85

library('vcd')

blood=read.table('mlbacDonor.k25.lz4.hist',skip=8)
#mda=read.table('mdaS23.k25.lz4.hist',skip=8)
#malbac=read.table('S01.k25.lz4.hist',skip=8)
indata=read.table(infile,skip=8)

getdata <- function(x,range) {
	thefreq <- x[,1]
	cnt <- x[,2]
	theratio <- x[,3]
	thedata <- cbind(thefreq,cnt)
	flag <- thedata[,1] <= range
	RVAL <- cbind(theratio[flag],thefreq[flag])
	RVAL
}

usedata1 <- getdata(blood,testXrange)
usedata2 <- getdata(indata,testXrange)
#usedata3 <- getdata(malbac,testXrange)

newfit <- function(x,y, method = c("ML", "MinChisq"), par = NULL) {
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
		newval <- as.vector(y[, 1])
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
	#expected <- sum(as.numeric(freq))/sum(p.hat) * p.hat
	expected <- p.hat
	expected <- newval
	df <- switch(method[1], MinChisq = {
		length(freq) + df
	}, ML = {
		sum(freq > 0) + df
	}, fixed = {
		c(length(freq), sum(freq > 0)) + df
	})
	RVAL <- list(observed = freq, count = count, fitted = expected,
		type = "poisson", method = method, df = df, par = par)
	class(RVAL) <- "goodfit"
	RVAL
}

fitres <- newfit(usedata2,usedata1, par = list(lambda = avgcvg))
#fitresMAL <- newfit(usedata3,usedata1, par = list(lambda = avgcvg))

sumres <- summary(fitres)
#sumresMAL <- summary(fitresMAL)

print(fitres)
print(sumres)

#./kmercmp.r mdaS23.k25.lz4.hist
#./kmercmp.r S01.k25.lz4.hist
