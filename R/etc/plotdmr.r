#!/usr/bin/env littler

da <- read.delim('dmrfile/Intron.dmr.out');
db <- read.delim('dmrfile/CDS.dmr.out');
dc <- read.delim('dmrfile/Upstream2k.dmr.out');
dd <- read.delim('dmrfile/CpGIsland.dmr.out');
de <- read.delim('dmrfile/Downstream2k.dmr.out');
df <- read.delim('dmrfile/5-UTR.dmr.out');
dg <- read.delim('dmrfile/3-UTR.dmr.out');

library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL;

#chr.exclude = c('chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15')

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
tracks.inside <- 7;
tracks.outside <- 0;

out.file <- "dmrCircosT.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
#png(file="dmrCircosT.png",2400,2400)

RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);
RCircos.Set.Plot.Area();

RCircos.Chromosome.Ideogram.Plot();
side <- "in";
data.col <- 6;
track.num <- 1;

cat("On Tumor",track.num,"\n");
RCircos.Histogram.Plot(da,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(db,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dc,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dd,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(de,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(df,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dg,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;

dev.off();

data.col <- 8;
track.num <- 1;
out.file <- "dmrCircosN.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
#png(file="dmrCircosN.png",2400,2400)

RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);
RCircos.Set.Plot.Area();

RCircos.Chromosome.Ideogram.Plot();

cat("On Normal",track.num,"\n");
RCircos.Histogram.Plot(da,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(db,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dc,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dd,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(de,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(df,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;
RCircos.Histogram.Plot(dg,data.col, track.num, side, min.value=0, max.value=1); track.num <- track.num+1;

dev.off();
