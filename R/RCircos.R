datN=read.delim('Normal.cg',header=F)
datB=read.delim('brother2.cg',header=F)

library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- NULL;
num.inside <- 5;
num.outside <- 0;
#chr.exclude = paste("chr", seq(5,22), sep="")

RCircos.Set.Core.Components(cyto.info,chr.exclude, num.inside, num.outside);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();

RCircos.Histogram.Plot(hist.data=datB,data.col=6, track.num=1, side="in");
RCircos.Histogram.Plot(hist.data=datN,data.col=6, track.num=2, side="in");

datN['D'] <- datB$V6-datN$V6
datN['X'] <- 20 * datN['D']
RCircos.Histogram.Plot(hist.data=datN,data.col=8, track.num=3, side="in");
RCircos.Histogram.Plot(hist.data=datN,data.col=7, track.num=5, side="in");

# https://github.com/cran/RCircos/blob/master/demo/RCircos.Demo.Human.R
data(RCircos.Histogram.Data);
RCircos.Histogram.Plot(hist.data=RCircos.Histogram.Data,data.col=4, track.num=6, side="in");
