library(tidyverse)
#library(Cairo)

if (0) {
	indat <- read.delim('e5.snp.txt', encoding="UTF-8")
	for (ncol in 1:nrow(indat)) {
		for (nrow in 5:ncol(indat)) {
			tmp <- unlist(strsplit(indat[ncol,nrow],','))
			indat[ncol,nrow] <- tmp[length(tmp)]
		}
	}
	write.csv(indat,'e5.depth.csv')
}

snpdat <- read.csv('e5.depth.csv', encoding="UTF-8")
plotdat <- pivot_longer(snpdat, cols = starts_with("rs"), names_to = "SNPs",values_to = "Depth")

PlotPointsize <- 12
cairo_pdf('boxType.pdf', pointsize = PlotPointsize)	# Not working on macOS, either cairo_pdf or CairoPDF.
ggplot(plotdat, aes(TypeCN, Depth)) + geom_boxplot(colour = rainbow(length(unique(plotdat$Type))))
dev.off()

pdf('boxSample.pdf', pointsize = PlotPointsize, width=32, height=16)
ggplot(plotdat, aes(ID, Depth)) + geom_boxplot(colour = rainbow(length(unique(plotdat$ID))))
dev.off()

pdf('boxTypeL10.pdf', pointsize = PlotPointsize)
ggplot(plotdat, aes(Type, Depth)) + geom_boxplot(colour = rainbow(length(unique(plotdat$Type)))) + scale_y_continuous(trans='log10')
dev.off()

pdf('boxSampleL10.pdf', pointsize = PlotPointsize, width=32, height=16)
ggplot(plotdat, aes(ID, Depth)) + geom_boxplot(colour = rainbow(length(unique(plotdat$ID)))) + scale_y_continuous(trans='log10')
dev.off()
