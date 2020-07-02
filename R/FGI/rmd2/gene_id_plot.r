#!/usr/bin/env littler

# Usage: Rscript gene_id_plot.R sample_report.csv

# 读入单个样本的报告表
#cmargs = commandArgs(trailingOnly = T)
sample_report_df = read.csv(argv[1])
#sample_report_df = read.csv("C:\\Users\\wangjingdong\\Desktop\\gene_id_rscript\\input_csv_report.csv")

# 准备数据框为目标画图形式
sample_report_df$End = sample_report_df$Pos
for (rsid in rownames(sample_report_df)){
  genotype = sample_report_df[rsid,"GT"]
  if (genotype == "Null"){sample_report_df[rsid, "Val"] = 0}
  else{
    allele_1 = substr(genotype, 1, 1)
    allele_2 = substr(genotype, 2, 2)
    if (allele_1 == allele_2){sample_report_df[rsid, "Val"] = 1}
    else{sample_report_df[rsid, "Val"] = 0.5}
  }
}

# 分成4个碱基文件
base_group = function(df, base){
  bool = c()
  for (i in 1:length(rownames(df))){
    geno = df[i,2]
    if (geno == "Null"){bool = append(bool,FALSE)}
    else{
      alleles = c(substr(geno, 1, 1), substr(geno, 2, 2))
      if (base %in% alleles){bool = append(bool,TRUE)}
      else{bool = append(bool,FALSE)}
    }
  }
  return(df[rownames(df)[bool],c("Chr","Pos","End","Val")])
}

matrix_A = base_group(sample_report_df, "A")
matrix_A$PlotColor = "#5050FF"

matrix_T = base_group(sample_report_df, "T")
matrix_T$PlotColor = "#CC9900"

matrix_G = base_group(sample_report_df, "G")
matrix_G$PlotColor = "#00C000"

matrix_C = base_group(sample_report_df, "C")
matrix_C$PlotColor = "#E00000"

#out.file <- "nCircos.pdf"
out.file <- argv[2]
#pdf(file=out.file, height=8, width=8)
png(file=out.file, height=2160, width=2160,pointsize=48)

library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram, tracks.inside=4, tracks.outside=0, chr.exclude=c("chrX", "chrY"))

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$track.background <- 'white'
rcircos.params$grid.line.color <- 'white'
rcircos.params$hist.width <- 160
rcircos.params$track.padding <- 0
RCircos.Reset.Plot.Parameters(rcircos.params)

RCircos.Set.Plot.Area(margins=0)
legend("center", legend=c("A", "C", "G", "T"), col=c("#5050FF", "#E00000","#00C000","#CC9900"), lwd=4)

RCircos.Chromosome.Ideogram.Plot()

RCircos.Histogram.Plot(hist.data=matrix_A, data.col=4,track.num=1, side="in",is.sorted=F)
RCircos.Histogram.Plot(hist.data=matrix_C, data.col=4,track.num=2, side="in",is.sorted=F)
RCircos.Histogram.Plot(hist.data=matrix_G, data.col=4,track.num=3, side="in",is.sorted=F)
RCircos.Histogram.Plot(hist.data=matrix_T, data.col=4,track.num=4, side="in",is.sorted=F)

dev.off()
