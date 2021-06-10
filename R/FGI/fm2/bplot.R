library(tidyverse)
snpdat <- read.csv('snps.csv')
plotdat <- pivot_longer(snpdat, cols = starts_with("DEP"), names_to = "SNPs",values_to = "Depth")

PlotPointsize <- 12
pdf('boxType.pdf', pointsize = PlotPointsize)
ggplot(plotdat, aes(Type, Depth)) + geom_boxplot()
dev.off()

pdf('boxSample.pdf', pointsize = PlotPointsize)
ggplot(plotdat, aes(ID, Depth)) + geom_boxplot()
dev.off()
