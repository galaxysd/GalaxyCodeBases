#!/usr/bin/env littler

snps=read.csv('outside.csv')
sdat=snps[,c(2,3,3,1,4,5,6,7)]
#         rsID   Chr       Pos      A      C      G      T    sum
#1  rs10500915 chr11  21753388 0.2891 0.4467 0.2643 0.0000 1.0001

out.file <- "nCircos.pdf";
pdf(file=out.file, height=8, width=8);

library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram, tracks.inside=4, tracks.outside=0, chr.exclude=c("chrX", "chrY")) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#RCircos.Tile.Plot(sdat, data.col=5, track.num=1, side="in",is.sorted=F)
# Scatter Heatmap Heatmap
RCircos.Heatmap.Plot(sdat, data.col=5, track.num=1, side="in",is.sorted=F)
RCircos.Heatmap.Plot(sdat, data.col=6, track.num=2, side="in",is.sorted=F)
RCircos.Heatmap.Plot(sdat, data.col=7, track.num=3, side="in",is.sorted=F)
RCircos.Heatmap.Plot(sdat, data.col=8, track.num=4, side="in",is.sorted=F)

dev.off();
