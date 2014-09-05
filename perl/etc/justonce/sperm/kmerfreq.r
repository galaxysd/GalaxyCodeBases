#!/usr/bin/env r
if (is.null(argv) | length(argv)<2) {
  cat("Usage: kmerfreq.r hist_file coverage\n")
  q()
}
infile <- argv[1]
avgcvg <- as.numeric(argv[2])

maxXrange <- 30
testXrange <- 20
#infile <- 'mlbacDonor.k25.lz4.hist'
#avgcvg <- 6.85

library('vcd')

aa=read.table(infile,skip=8)
thefreq <- aa[,1]
cnt <- aa[,2]
theratio <- aa[,3]

thedata <- cbind(thefreq,cnt)
flag <- thedata[,1] <= testXrange
#usedata <- cbind(cnt[flag],thefreq[flag])
usedata <- cbind(theratio[flag],thefreq[flag])

#fitres <- goodfit(usedata,type='poisson', par = list(lambda = avgcvg))

#sumres <- summary(fitres)

newfit <- function(x, method = c("ML", "MinChisq"), par = NULL) {
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
	#expected <- sum(as.numeric(freq))/sum(p.hat) * p.hat
	expected <- p.hat
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

fitres <- newfit(usedata, par = list(lambda = avgcvg))

sumres <- summary(fitres)

print(fitres)
print(sumres)


tiff(filename = paste0(argv[1],".tiff"), compression="lzw", width = 683, height = 683)#, units = "px", pointsize = 52)
par(mar=c(5, 8, 2, 0.5),ps=20,family='sans')

yy=aa[,3]
xx=aa[,1]
#xx <- usedata[,2]	# count
#yy <- usedata[,1]	# freq

mean2=avgcvg
zz=dpois(xx,mean2);
#lines(xx,zz,xlim=c(0,60),type='l');

#barplot(yy,xlim=c(0,20),ylim=c(0,0.212),xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
plot(xx,yy,type='h',lwd=20,xlim=c(1,maxXrange),ylim=c(0,0.212),
	main = paste0("Pearson goodness of fit, p = ",sumres[1,3]),
	xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
# dpois(3,3.62)=0.2117524
legend("topright",pch=c(15,-1),lty=c(-1,1),col=c("navy","red"), x.intersp = 1, y.intersp = 2,cex=1,lwd=4,
	legend= c(argv[1],paste0("lamda=",argv[2])))
axis(1, at = seq(0, maxXrange, by = 5),lwd=3,cex=1)
lines(xx,zz,xlim=c(0,maxXrange),type='l',col='red',lwd=2)
#lines(fitres$count,fitres$observed,xlim=c(0,60),type='l',col='blue')
lines(fitres$count,fitres$fitted,xlim=c(0,60),type='l',col='red',lwd=6)
dev.off()


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
