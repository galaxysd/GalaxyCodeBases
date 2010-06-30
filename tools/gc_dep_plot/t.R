a=read.delim("./HUMgscRFDDIAAPEI-12.soaplst.dat");
pdf("./HUMgscRFDDIAAPEI-12.pdf",width=8,height=6,paper="special");
b=boxplot(a,plot=F);
m=ceiling(max(b$stats,na.rm=T));
b$stats;

m=3;
bxp(b,axes=F,outline=F,xlab='GC%',ylab='Depth',ylim=c(0,m),xlim=c(0,90)); 
axis(1,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=c(0,10,20,30,40,50,60,70,80,90,100));
axis(2,las=2);
dev.off();
