R CMD INSTALL --no-docs /ifs2/POPULATION/PRJECT/Rice_QQ/sprice/genotyping/qtl/qtl_1.18-7.tar.gz

library('qtl')
dat=read.cross('csvsr','.','gen_rot.csv','phe_rot.csv')

# 查看表型分布
par(mfrow=c(2,1))
for(i in 1:2)
  plot.pheno(dat,i)

# 查看极端表型值
pairs(jitter( as.matrix(dat$pheno[,1:2]) ), cex=0.6, las=1)

# 查看表型与个体编号的相关性
par(mfrow=c(1,2), las=1, cex=0.8)
means <- apply(dat$pheno[,1:2], 1, mean)
plot(means)
plot(sample(means), xlab="Random index", ylab="means")

# 查看个体间表型相似度分布
cg <- comparegeno(dat)
par(mar=c(4.1,1.1,0.1,1.1))
hist(cg, breaks=seq(0, 1, len=201), main="", yaxt="n", ylab="",
     xlab="Proportion of identical genotypes")
rug(cg)
which(cg > 0.9, arr.ind=TRUE)

# 计算个体交换次数
nxo <- countXO(dat)
plot(nxo, ylab="No. crossovers")
# nxo[nxo>1000]

z <- plot.info(dat, step=0)
z[ z[,1]==2, ]


out <- scanone(dat, pheno.col=2)
write.table(out,'height.txt')
png('height.png',16000,1200)
plot(out, ylab="LOD score",main='Raw SNP Markers')
dev.off()

dat0 <- dat
dat <- calc.genoprob(dat, step=1, error.prob=0.01)
out.np <- scanone(dat, model="np", pheno.col=2)
write.table(out.np,'height.np.txt')
png('height.np.png',16000,1200)
plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
dev.off()


#operm.np <- scanone(listeria, model="np", n.perm=1000, perm.Xsp=TRUE)
#operm.np <- scanone(dat, model="np", n.perm=100, pheno.col=2)
#summary(out.np, perms=operm.np, alpha=0.05, pvalues=TRUE)
#write.table(operm.np,'height.operm.np.txt')

out.mr <- scanone(dat, method="mr")
out.2p <- scanone(dat, pheno.col=1, method="mr")





library('qtl')
file='av1'
dat=read.cross('csvsr','.',paste(sep='',file,'_rot.csv'),'2nd_phe.csv')

out <- scanone(dat, pheno.col=1)
write.table(out,paste(sep='',file,'_chlorophyl_content.txt'))
png(paste(sep='',file,'_chlorophyl_content.png'),16000,1200)
plot(out, ylab="LOD score",main='Raw SNP Markers')
dev.off()

dat0 <- dat
dat <- calc.genoprob(dat, step=1, error.prob=0.01)
out.np <- scanone(dat, model="np", pheno.col=1)
write.table(out.np,paste(sep='',file,'_chlorophyl_content.np.txt'))
png(paste(sep='',file,'_chlorophyl_content.np.png'),16000,1200)
plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
dev.off()

out <- scanone(dat0, pheno.col=2)
write.table(out,paste(sep='',file,'_heading_days.txt'))
png(paste(sep='',file,'_heading_days.png'),16000,1200)
plot(out, ylab="LOD score",main='Raw SNP Markers')
dev.off()
out.np <- scanone(dat, model="np", pheno.col=2)
write.table(out.np,paste(sep='',file,'_heading_days.np.txt'))
png(paste(sep='',file,'_heading_days.np.png'),16000,1200)
plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
dev.off()

out <- scanone(dat0, pheno.col=3)
write.table(out,paste(sep='',file,'_plant_height2.txt'))
png(paste(sep='',file,'_plant_height2.png'),16000,1200)
plot(out, ylab="LOD score",main='Raw SNP Markers')
dev.off()
out.np <- scanone(dat, model="np", pheno.col=3)
write.table(out.np,paste(sep='',file,'_plant_height2.np.txt'))
png(paste(sep='',file,'_plant_height2.np.png'),16000,1200)
plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
dev.off()

out <- scanone(dat0, pheno.col=4)
write.table(out,paste(sep='',file,'_tiller_number2.txt'))
png(paste(sep='',file,'_tiller_number2.png'),16000,1200)
plot(out, ylab="LOD score",main='Raw SNP Markers')
dev.off()
out.np <- scanone(dat, model="np", pheno.col=4)
write.table(out.np,paste(sep='',file,'_tiller_number2.np.txt'))
png(paste(sep='',file,'_tiller_number2.np.png'),16000,1200)
plot(out.np, ylab="LOD score", alternate.chrid=TRUE,main='HMM filtered')
dev.off()



./qtl.R <gen_rot.csv> <phe_rot.csv> <output_prefix>
v1_rot.csv
av1_rot.csv

fv1_phe.csv
2nd_phe.csv
3rd_phe.csv
