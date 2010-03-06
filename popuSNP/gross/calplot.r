a=read.delim('./outpoly/mix.allpoly',header=F);
p=0.05;

theX=qnorm(p*2);

PiR=a$V2/a$V5;
TDD=a$V4;
PiRi=PiR[!is.infinite(PiR)];
TDDi=TDD[!is.infinite(PiR)];	# not necessary
mp=mean(PiRi,na.rm=T);
sdp=sd(PiRi,na.rm=T);
lowp=mp+theX*sdp;

mt=mean(TDDi,na.rm=T);
sdt=sd(TDDi,na.rm=T);
lowt=mt+theX*sdt;


PiRH=hist(PiR,br='scott',plot=FALSE);
TDDH=hist(TDD,br='scott',plot=FALSE);
top1=max(PiRH$counts)/sum(PiRH$counts);
top2=max(TDDH$counts)/sum(TDDH$counts);

def.par <- par(no.readonly = TRUE)
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1),c(1,3),TRUE);
par(mar=c(5,5,0,0));
plot(PiR,TDD,xlim=c(0,5),ylim=c(-3,4),pch="+",xlab=expression(theta[pi][","][cultivate]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[cultivate]));

#plot(function(x)dnorm(x), -6, 6)

x1=c(lowp,lowp);
y1=c(-3.3,4.3);
lines(x1,y1,col="RED");
x1=c(-0.3,5.3);
y1=c(lowt,lowt);
lines(x1,y1,col="RED");
par(mar=c(0,5,5,1));
barplot(PiRH$counts/sum(PiRH$counts),axes=FALSE, ylim=c(0,top1), space=0, col="BLUE", ylab="Proportion",cex.lab=1.5);
ypos = seq(0,0.1,by=0.05)
axis(side=2, ypos, labels=FALSE)
mtext(ypos, side=2,at=ypos, line=1, cex=1.2)
threshold_x2=4.31;
x2=c(threshold_x2,threshold_x2);
y2=c(0,top1);
lines(x2,y2,col="RED");
par(mar=c(5,0,1,3));
barplot(TDDH$counts/sum(TDDH$counts),axes=FALSE,xlim=c(0,top2),space=0,horiz=TRUE,col="GREEN",xlab="Proportion",cex.lab=1.5);
ypos = seq(0,0.03,by=0.015)
axis(side=1, ypos, labels=FALSE)
mtext(ypos, side=1,at=ypos,line=1, cex=1.2)
threshold_y3=34.09;
x3=c(0,top2);
y3=c(threshold_y3,threshold_y3);
lines(x3,y3,col="RED");

dev.off();
par(def.par)

myhist = function(data)
{
	postscript("tempNew.ps");
	opar = par()
	a = data
	PiR=a$V2/a$V5
	TDD=a$V4;
	x_step=0.05;
	y_step=0.05;
	PiRH=hist(PiR,br=seq(0,5,x_step),plot=FALSE);
	TDDH=hist(TDD,br=seq(-3,4,y_step),plot=FALSE);
	top1=max(PiRH$counts/sum(PiRH$counts));
	top2=max(TDDH$counts/sum(TDDH$counts));
	layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1),c(1,3),TRUE);
	par(mar=c(5,5,0,0));
	plot(PiR,TDD,xlim=c(0,5),ylim=c(-3,4),pch="+",xlab=expression(theta[pi][","][domesticated]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[domesticated]));
	lowp=0.2;
	lowt=-1.37;
	x1=c(lowp,lowp);
	y1=c(-3.3,4.3);
	lines(x1,y1,col="RED");
	x1=c(-0.3,5.3);
	y1=c(lowt,lowt);
	lines(x1,y1,col="RED");
	par(mar=c(0,5,5,1));
	barplot(PiRH$counts/sum(PiRH$counts),axes=FALSE, ylim=c(0,top1), space=0, col="BLUE", ylab="Proportion",cex.lab=1.5);
	ypos = seq(0,0.1,by=0.05)
	axis(side=2, ypos, labels=FALSE)
	mtext(ypos, side=2,at=ypos, line=1, cex=1.2)
	threshold_x2=4.31;
	x2=c(threshold_x2,threshold_x2);
	y2=c(0,top1);
	lines(x2,y2,col="RED");
	par(mar=c(5,0,1,3));
	barplot(TDDH$counts/sum(TDDH$counts),axes=FALSE,xlim=c(0,top2),space=0,horiz=TRUE,col="GREEN",xlab="Proportion",cex.lab=1.5);
	ypos = seq(0,0.03,by=0.015)
	axis(side=1, ypos, labels=FALSE)
	mtext(ypos, side=1,at=ypos,line=1, cex=1.2)
	threshold_y3=34.09;
	x3=c(0,top2);
	y3=c(threshold_y3,threshold_y3);
	lines(x3,y3,col="RED");
	dev.off();
	#ypos = seq(8,18,by=2)
	#par(mar=c(4.5, 4.5, 2.5, 2.5))
	#plot(rt, type='p',col='blue' ,pch=19, cex=0.6, xaxt="n",yaxt="n", xlab="", ylab="",ylim=c(8,18) )
	#axis(side=1, xpos, tcl=0.2, labels=FALSE)
	#axis(side=2, ypos, tcl=0.2, labels=FALSE)
	#mtext("Samples",side=1, line=2, at=median(xpos), cex=1.9 )
	#mtext("Depth of target region",side=2, line=3, at=median(ypos), cex=1.8 )
	#mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
	#mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	#par(opar)
	#dev.off()
}
