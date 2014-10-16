D__myfit <- 1

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
			p.hat <- diff(c(0, ppois(count[-n], lambda = x),1)) # DO NOT USE THIS NOW.
			expected <- sum(freq) * p.hat
			sum((freq - expected)^2/expected)
		}
		par <- optimize(chi2, range(count))$minimum
	}
	par <- list(lambda = par)
	p.hat <- dpois(count, lambda = par$lambda)

	zoomratio <- sum(freq) / sum(p.hat)
	expected <- p.hat * zoomratio
	print(c(sum(freq),sum(p.hat),zoomratio))
	#print(cbind(p.hat,expected))

	df <- switch(method[1], MinChisq = {
		length(freq) + df
	}, ML = {
		sum(freq > 0) + df
	}, fixed = {
		c(length(freq), sum(freq > 0)) + df
	})
	RVAL <- list(observed = freq, count = count, fitted = expected,
		type = "poisson", method = method, df = df, par = par, zoom = zoomratio)
	class(RVAL) <- "myfit"
	RVAL
}



scale2pois <- function(x, lambda, cntmin = NULL, cntmax = NULL) {
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

	p.hat <- dpois(count, lambda = lambda)
	if (is.null(cntmin) | is.na(cntmin))
		cntmin <- min(count)
	if (is.null(cntmax) | is.na(cntmax))
		cntmax <- max(count)
	flag <- (count >= cntmin) & (count <= cntmax)
	#print(c(cntmin,cntmax,count,flag))
	zoomratio <- sum(p.hat[flag]) / sum(freq[flag])
	scaled <- freq * zoomratio
	#expected <- p.hat / zoomratio
	#print(c(sum(freq),sum(scaled),sum(p.hat),sum(expected),zoomratio))
	#print(cbind(freq,scaled, p.hat,expected))

	RVAL <- list(ret = cbind(scaled,count), zoom = zoomratio, cntmin = cntmin, cntmax = cntmax)
	RVAL
}

cdf2pdf <- function(x, rescale = FALSE) {
	mylen <- length(x) - 1
	#RVAL <- diff(c(0, x, 1))[1:mylen]
	RVAL <- diff(c(0, x), 1)[1:mylen]
	if (rescale) {
		s <- sum(RVAL)
		RVAL <- RVAL/s
	}
	RVAL
}

# http://www.csee.ogi.edu/~gormanky/papers/handbook/ky.R
pdf2cdf <- function(x, complete = FALSE, rescale = FALSE) {
    # make a sorted PDF into a sorted CDF
	mylen <- length(x)
	RVAL <- x
	if (complete | rescale) {
		RVAL <- c(x,x[mylen])
		mylen <- 1 + mylen
		#length(RVAL) <- mylen
	}
	for (i in 2:mylen) {
		RVAL[i] <- RVAL[i] + RVAL[i-1]
	}
	if (rescale) {
		RVAL <- RVAL/RVAL[mylen]
	}
	return(RVAL)
}
