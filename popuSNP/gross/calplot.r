#a=read.delim('dat.tsv',header=F);
a=read.delim('./outpoly/mix.allpoly',header=F);
#a=read.delim('E:\\BGI\\toGit\\popuSNP\\gross\\chr3.polymorphism',header=F);
#a=read.delim('E:\\BGI\\toGit\\popuSNP\\gross\\dat10k.tsv',header=F);
p=0.01;

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

cat(lowp,lowt);

PiRH=hist(PiR,br='scott',plot=FALSE);
TDDH=hist(TDD,br='scott',plot=FALSE);
top1=ceiling(max(PiRH$counts)/sum(PiRH$counts)/0.02)*2;
top2=ceiling(max(TDDH$counts)/sum(TDDH$counts)/0.02)*2;

png('out.png', 1000, 1000, pointsize=12,res=96);
#def.par <- par(no.readonly = TRUE) # save default, for resetting...
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1),c(1,3),TRUE);
par(mar=c(5,5,0,0));
plot(PiR,TDD,xlim=c(0,5),ylim=c(-3,4),pch="+",xlab=expression(theta[pi][","][cultivate]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[cultivate]));
#plot(PiR,TDD,ylim=c(-3,4),log='x',pch="+",xlab=expression(theta[pi][","][cultivate]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[cultivate]));

#plot(function(x)dnorm(x), -6, 6)

x1=c(lowp,lowp);
y1=c(-9993.3,9994.3);
lines(x1,y1,col="RED");
x1=c(-9990.3,9995.3);
y1=c(lowt,lowt);
lines(x1,y1,col="RED");

par(mar=c(0,5,5,0));
Pb=PiRH$breaks[2]-PiRH$breaks[1];
barplot(100*PiRH$counts/sum(PiRH$counts),width=Pb,xlim=c(0,5),ylim=c(0,top1), space=0, col="BLUE", ylab="Proportion (%)",cex.lab=1.1,cex.axis=1);
#threshold_x2=4.31;
x2=c(lowp,lowp);
y2=c(0,top1);
lines(x2,y2,col="RED");

par(mar=c(5,0,0,3));
Tb=TDDH$breaks[2]-TDDH$breaks[1];
barplot(100*TDDH$counts/sum(TDDH$counts),width=Tb,ylim=c(-3,4)+3,xlim=c(0,top2),space=0,horiz=TRUE,col="GREEN",xlab="Proportion (%)",cex.lab=1.1,cex.axis=1);
#threshold_y3=34.09;
x3=c(0,top2);
y3=c(lowt,lowt)+3;
lines(x3,y3,col="RED");

dev.off();
#par(def.par)
