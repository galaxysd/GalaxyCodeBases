#a=read.delim('dat.tsv',header=F);
a=read.delim('./outpoly/mix.mpoly',header=F);
#a=read.delim('E:\\BGI\\toGit\\popuSNP\\gross\\chr3.polymorphism',header=F);
#a=read.delim('E:\\BGI\\toGit\\popuSNP\\gross\\dat10k.tsv',header=F);
a=na.omit(a);
p=0.01;

theX=qnorm(p*2);

PiR=a$V2/a$V5;
TDD=a$V4;
PiRi=PiR[!is.infinite(PiR)];
TDDi=TDD[!is.infinite(PiR)];	# not necessary
pp=quantile(PiRi,c(0.01,1));
dd=quantile(TDDi,c(0.01,1));
ppp=quantile(PiRi,c(0,.9999));
#ddd=quantile(TDDi,c(.005,.995));
ppp=c(0,ceiling(ppp[2]));
#ddd=c(floor(ddd[1]),ceiling(ddd[2]));
ddd=c(floor(min(TDDi)),ceiling(max(TDDi)));

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
plot(PiR,TDD,xlim=ppp,ylim=ddd,pch="+",xlab=expression(theta[pi][","][cultivate]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[cultivate]));
#plot(PiR,TDD,ylim=c(-3,4),log='x',pch="+",xlab=expression(theta[pi][","][cultivate]/theta[pi][","][wild]),cex.lab=1.5,cex.axis=1.5,ylab=expression(Tajima ~~ italic(D)[cultivate]));

#plot(function(x)dnorm(x), -6, 6)

x1=c(lowp,lowp);
x1=c(pp[1],pp[1]);

y1=c(-9993.3,9994.3);
lines(x1,y1,col="RED");
x1=c(-9990.3,9995.3);
y1=c(lowt,lowt);
y1=c(dd[1],dd[1]);

lines(x1,y1,col="RED");

par(mar=c(0,5,5,0));
Pb=PiRH$breaks[2]-PiRH$breaks[1];
barplot(100*PiRH$counts/sum(PiRH$counts),width=Pb,xlim=ppp,ylim=c(0,top1), space=0, col="BLUE", ylab="Proportion (%)",cex.lab=1.1,cex.axis=1);
#threshold_x2=4.31;
x2=c(lowp,lowp);
x2=c(pp[1],pp[1]);

y2=c(0,top1);
lines(x2,y2,col="RED");

par(mar=c(5,0,0,3));
Tb=TDDH$breaks[2]-TDDH$breaks[1];
barplot(100*TDDH$counts/sum(TDDH$counts),width=Tb,ylim=ddd-ddd[1],xlim=c(0,top2),space=0,horiz=TRUE,col="GREEN",xlab="Proportion (%)",cex.lab=1.1,cex.axis=1);
#threshold_y3=34.09;
x3=c(0,top2);
y3=c(lowt,lowt)-ddd[1];
y3=c(dd[1],dd[1])-ddd[1];

lines(x3,y3,col="RED");

dev.off();
#par(def.par)

#perl -lane '@a=split /\t/;next if $a[4]==0;$p=$a[1]/$a[4];$d=$a[3];print if $p<=0.02321981 and $d<=-1.995619;' mix_5k_cn116_f01.mpoly > mix_5k_cn116_f01.filtered
