#!/usr/bin/env littler

########################################################################################################
#
# Usage: Rscript gene_id_run.R -i input_samples.csv -f allele_freq.csv -d database.csv -m input_info.txt
#
########################################################################################################

# -i files/msSNP.csv -m files/meta.txt -d files/database.csv -f files/allele_freq.csv -p out

# 设置数值显示位数
options(scipen = 200)

# 读取命令行参数
library(getopt)
spec = matrix(
  c("input","i",2,"character","the input genotype csv of suspect",
    "allele_freq","f",2,"character","allele frequence csv of SNPs",
    "database","d",2,"character","the database csv stored SNPs genotype",
	"outpath","p",2,"character","out dir.",
    "input_info","m",2,"character","meta infomation of input samples"),
  byrow = TRUE, ncol = 5
)

opt = getopt(spec = spec)

library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)

# 定义一个计算两个样本相似度的核心函数
calculate = function(input_geno_df, sample_geno_df, allele_freq_df){
  likehood = 1
  for (rsid in colnames(input_geno_df)){
    ig = input_geno_df[1,rsid]
    sg = sample_geno_df[1,rsid]
    if ((ig == 'Null') || (sg == 'Null')){pm_value = 1}
    else{
      if (ig != sg){pm_value = 0.0001}
      if (ig == sg){
        iw = c(substr(ig,1,1), substr(ig,2,2))
        sw = c(substr(sg,1,1), substr(sg,2,2))
        if (iw[1] == iw[2]){pm_value = 1 / ((allele_freq_df[rsid,iw[1]])**2)}
        if (iw[1] != iw[2]){
          if ((allele_freq_df[rsid,iw[1]]==0.01) || (allele_freq_df[rsid,iw[2]]==0.01)){pm_value = 100}
          else{pm_value = 1 / (2 * allele_freq_df[rsid,iw[1]] * allele_freq_df[rsid,iw[2]])}
        }
      }
    }
    sample_geno_df[1,rsid] = paste(sg, pm_value)
    likehood = likehood * pm_value
  }
  sample_geno_df$tpm = likehood
  return(sample_geno_df)
}

# 读入文件
input_geno_df = read.csv(opt$input, row.names = 1)
#allele_freq_df = read.csv(opt$allele_freq, row.names = 1)
#database_df = read.csv(opt$database, row.names = 1)
input_info_df = suppressWarnings(read.delim(opt$input_info, fileEncoding = "UTF-16", row.names = 1))
#input_geno_df = read.csv("C:\\Users\\wangjingdong\\Desktop\\gene_id_rscript\\gene_id_rscript\\data\\input_csv.csv", row.names = 1)
allele_freq_df = read.csv("files/allele_freq.csv", row.names = 1)
database_df = read.csv("files/database.csv", row.names = 1)
#input_info_df = read.delim("C:\\Users\\wangjingdong\\Desktop\\gene_id_rscript\\gene_id_rscript\\data\\input_info.txt", fileEncoding = "UTF-16", row.names = 1)

# 创建样本比对结果输出文件夹
input_file_name = strsplit(basename(opt$input),".", fixed = T)[[1]][1]
#input_file_name = strsplit(basename("C:\\Users\\wangjingdong\\Desktop\\gene_id_rscript\\gene_id_rscript\\data\\input_csv.csv"),".", fixed = T)[[1]][1]
#out_dir_name = paste(input_file_name, "_out", sep = '')
#out_dir_name = paste(opt$outpath)
out_dir_name <- 'out'
dir.create(out_dir_name,showWarnings=F)
dir.create('pdf',showWarnings=F)

# 对数据预处理，将等位基因频率0替换成0.01，将纯合‘A’替换成’AA‘
allele_freq_df[allele_freq_df==0] = 0.01
input_geno_df[input_geno_df == "TRUE"] = "T"
for (r in rownames(input_geno_df)){
  for (c in colnames(input_geno_df)){
    if (nchar(input_geno_df[r,c]) == 1){
      input_geno_df[r,c] = paste(input_geno_df[r,c], input_geno_df[r,c], sep = '')
    }
  }
}

for (r in rownames(database_df)){
  for (c in colnames(database_df)){
    if (nchar(database_df[r,c]) == 1){
      database_df[r,c] = paste(database_df[r,c], database_df[r,c], sep = '')
    }
  }
}

# 计算每个样本与数据库的比对结果，并返回比对总表、基因型报告两个表格
NeedRerun = c()
for (suspect in rownames(input_geno_df)){
  match_file = paste(suspect, "_match.csv", sep = '')
  report_file = paste(suspect, "_report.csv", sep = '')
  suspect_geno_df = input_geno_df[suspect,]
  cat('Processing',suspect,':')
  # 计算该样本与数据库中每个样本的比对结果
  samples = rownames(database_df)
  sample_geno_df_one = calculate(suspect_geno_df, database_df[samples[1],], allele_freq_df)
  for (sample_name in samples[c(2:length(samples))]){
    sample_geno_df_two = calculate(suspect_geno_df, database_df[sample_name,], allele_freq_df)
    sample_geno_df_one = rbind(sample_geno_df_one, sample_geno_df_two)
  }
  sorted_all_result_df = sample_geno_df_one[order(sample_geno_df_one[,"tpm"], decreasing = T), ]
  # 对匹配结果进行分级，做判断陈述
  QCstate = 'PASS'
  if (sorted_all_result_df[1,"tpm"] < 0.0001){statement = "无匹配对象"}else{
    if (sorted_all_result_df[1,"tpm"] <= 8000000000){
      mid_match = subset(sorted_all_result_df, (tpm <= 8000000000) && (tpm >= 0.0001))
      samples = paste(rownames(mid_match), collapse = ',')
      statement = paste("与样本", samples, "可能一致, 应增加其它检测")
      QCstate = 'FAIL'
      NeedRerun <- c(NeedRerun,suspect)
    }else{
      top_match = subset(sorted_all_result_df, tpm > 8000000000)
      samples = paste(rownames(top_match), collapse = ',')
      statement = paste("与样本", samples, "一致")
    }
  }
  sorted_all_result_df$tpm = format(sorted_all_result_df$tpm, scientific = T)
  # 输出单个样本比对总表，要根据系统调整文件名分隔符，‘/’或者‘\\’
  write.csv(sorted_all_result_df, file.path(out_dir_name, match_file), fileEncoding = "UTF-8")

  # 输出单个样本的基因分型报告
  best_rate_person = rownames(sorted_all_result_df)[1]
  best_rate = sorted_all_result_df[best_rate_person, "tpm"]

  # 合成初步的报告表格
  trans_suspect_geno_df = as.data.frame(t(suspect_geno_df))
  names(trans_suspect_geno_df) = "GT"
  report_df = cbind(trans_suspect_geno_df, allele_freq_df[,c("Number","Chr","Pos","A","C","G","T")])

  # 计算人群基因型频率百分比
  for (rsid in rownames(report_df)){
    genotype = report_df[rsid,"GT"]
    if (genotype == "Null"){report_df[rsid, "GT%"] = "0%"}
    else{
      allele_1 = substr(genotype, 1, 1)
      allele_2 = substr(genotype, 2, 2)
      if (allele_1 == allele_2){report_df[rsid, "GT%"] = paste(format((report_df[rsid, allele_1]) ** 2 * 100, digits = 4),"%", sep = '')}
      else{report_df[rsid, "GT%"] = paste(format(2 * report_df[rsid, allele_1] * report_df[rsid, allele_2] * 100, digits = 4), "%", sep = '')}
    }
  }

  report_df$rsID = rownames(report_df)
  report_df[rownames(report_df)[1],"report"] = suspect
  report_df[rownames(report_df)[2],"report"] = best_rate
  report_df[rownames(report_df)[3],"report"] = best_rate_person
  report_df[rownames(report_df)[4],"report"] = input_info_df[suspect, "Name"]
  report_df[rownames(report_df)[5],"report"] = input_info_df[suspect, "Sex"]
  report_df[rownames(report_df)[6],"report"] = input_info_df[suspect, "Birth"]
  report_df[rownames(report_df)[7],"report"] = input_info_df[suspect, "Date"]
  report_df[rownames(report_df)[8],"report"] = input_info_df[suspect, "Tel"]
  report_df[rownames(report_df)[9],"report"] = input_info_df[suspect, "Address"]
  report_df[rownames(report_df)[10],"report"] = statement
  report_df[rownames(report_df)[11],"report"] = QCstate
  output_report_df = report_df[,c("Number","GT","GT%","rsID","Chr","Pos","report")]

  # 输出单个样本基因型报告表格，要根据系统调整文件名分隔符，‘/’或者‘\\’
  write.csv(output_report_df, file.path(out_dir_name, report_file), row.names = F, fileEncoding = "UTF-8")

  # Plot
  sample_report_df <- output_report_df
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
  png_file <- paste(suspect, "_report.png", sep = '')
  png(file=file.path(out_dir_name, png_file), height=2160, width=2160,pointsize=48)

  #RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram, tracks.inside=4, tracks.outside=0, chr.exclude=c("chrX", "chrY"))
  suppressMessages(
	  RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram, tracks.inside=4, tracks.outside=0, chr.exclude=c("chrX", "chrY"))
  , classes = "message")

  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$track.background <- 'white'
  rcircos.params$grid.line.color <- 'white'
  rcircos.params$hist.width <- 160
  rcircos.params$track.padding <- 0
  suppressMessages(
	  RCircos.Reset.Plot.Parameters(rcircos.params)
  , classes = "message")

  RCircos.Set.Plot.Area(margins=0)
  legend("center", legend=c("A", "C", "G", "T"), col=c("#5050FF", "#E00000","#00C000","#CC9900"), lwd=16)

  RCircos.Chromosome.Ideogram.Plot()

  RCircos.Histogram.Plot(hist.data=matrix_A, data.col=4,track.num=1, side="in",is.sorted=F)
  RCircos.Histogram.Plot(hist.data=matrix_C, data.col=4,track.num=2, side="in",is.sorted=F)
  RCircos.Histogram.Plot(hist.data=matrix_G, data.col=4,track.num=3, side="in",is.sorted=F)
  RCircos.Histogram.Plot(hist.data=matrix_T, data.col=4,track.num=4, side="in",is.sorted=F)

  dev.off()

  png_trimed <- paste(suspect, ".png", sep = '')
  imagemagickCLI <- paste('convert -trim +repage -bordercolor none -border 0x20',file.path(out_dir_name, png_file),file.path(out_dir_name, png_trimed))
  system(imagemagickCLI)

  thePrefix <- paste(suspect, "_report", sep = '')
  invisible(capture.output(
	  rmarkdown::render('20200608rep.Rmd',output_dir=out_dir_name,output_file=paste0(thePrefix,'.pdf'),quiet=TRUE,clean=TRUE)
  ))

  mergeCLI <- paste('qpdf --empty --pages inc/Cover.pdf ',file.path(out_dir_name, paste0(thePrefix,'.pdf')),' -- ', file.path('pdf', paste0(thePrefix,'.pdf')) )
  system(mergeCLI)

  cat(' done.\n')
}

qcfilecon <- file(file.path(out_dir_name, '_QC.txt'),'wt',encoding='UTF-8')
if (length(NeedRerun)) {
	QCstate <- paste(sep="\n",'下列样品需要重测：',NeedRerun,'')
} else {
	QCstate <- '这批样品正常。\n'
}
cat('\n=================\n',QCstate,'=================\n',sep='')
writeLines(QCstate,con=qcfilecon)
close(qcfilecon)
