#!/usr/bin/env r
if (is.null(argv) | length(argv)<2) {
  cat("Usage: kmerfreq.r hist_file coverage\n")
  q()
}
testXrange <- as.numeric(argv[1])
infile <- argv[2]
avgcvg <- as.numeric(argv[3])

maxXrange <- 25
#testXrange <- 20
#infile <- 'blood.fq.k25.gz.hist'
#avgcvg <- 3.27885

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
	#expected <- sum(freq)/sum(p.hat) * p.hat
	#expected <- p.hat
	expected <- p.hat * sum(freq) / sum(p.hat)
	print(c(sum(freq),sum(p.hat),sum(freq) / sum(p.hat)))
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
		type = "poisson", method = method, df = df, par = par)
	class(RVAL) <- "goodfit"
	RVAL
}

selfres <- newfit(usedata, par = list(lambda = avgcvg))
selfsum <- summary(selfres)
selfchi <- chisq.test(cbind(selfres$observed,selfres$fitted))
print(selfres)
print(selfsum)
print(selfchi)

newres<-list(observed = fitres$observed[c(-1,-2)], count = fitres$count[c(-1,-2)], fitted = fitres$fitted[c(-1,-2)],
		type = "poisson", method = fitres$method, df = fitres$df - 2, par = fitres$par)
class(newres) <- "goodfit"

sumres <- summary(newres)

print(newres)
print(sumres)

print(cbind(newres$observed,newres$fitted,(newres$observed-newres$fitted)/newres$fitted, selfres$fitted,(selfres$observed-selfres$fitted)/selfres$fitted))

tiff(filename = paste0(infile,".tiff"), compression="lzw", width = 683, height = 683)#, units = "px", pointsize = 52)
par(mar=c(5, 8, 2, 0.5),ps=20,family='sans')

yy=aa[,3]
xx=aa[,1]
#xx <- usedata[,2]	# count
#yy <- usedata[,1]	# freq

mean2=avgcvg
xxx=seq(0,max(xx))
zz=dpois(xxx,mean2);
maxY=max(zz,yy)
maxY<-0.3
#lines(xx,zz,xlim=c(0,60),type='l');

#barplot(yy,xlim=c(0,20),ylim=c(0,0.212),xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
plot(xx,yy,type='h',lwd=20,xlim=c(1,maxXrange),ylim=c(0,maxY),
	main = paste0("Pearson goodness of fit, p=",sumres[1,3]),
	xlab="K-mer count",ylab='Ratio',col='navy',cex.lab=1)
# dpois(3,3.62)=0.2117524
legend("topright",pch=c(15,-1),lty=c(-1,1),col=c("navy","red"), x.intersp = 1, y.intersp = 2,cex=1,lwd=4,
	legend= c(argv[1],paste0("lamda=",argv[2])))
axis(1, at = seq(0, maxXrange, by = 5),lwd=3,cex=1)
#lines(fitres$count,fitres$observed,xlim=c(0,60),type='l',col='blue')
lines(fitres$count,fitres$fitted,xlim=c(0,60),type='l',col='red',lwd=6)
lines(selfres$count,selfres$fitted,xlim=c(0,60),type='l',col='green',lwd=6)
lines(xxx,zz,xlim=c(0,maxXrange),type='l',col='black',lwd=2)
dev.off()

print(selfchi)

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
