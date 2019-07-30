library(ggplot2)
library(splines)
args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]

#pdf(file=output,height=6.18,width=10)
pdf(output)
raw.tab <- read.table(file=input,sep ="\t",header = T)
input.tab <- subset(raw.tab, select = -1 )
name <- names(input.tab)


for (i in 1:ncol(input.tab)){
	print(
		ggplot(data=input.tab,aes(x=input.tab[,i])) +
		geom_histogram(position="dodge", bins = 200,na.rm = TRUE, color="darkblue", fill="lightblue") +
		labs(title=name[i],x="cffDNA(%)",y="Number of SNP",size=20) +
#		theme_bw() +
#		theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
		scale_y_continuous(breaks=seq(0,10000,10)) +
		scale_x_continuous(limits=c(0,50),breaks=seq(0,100,10)) +
		theme(axis.title.x = element_text(size = 10, color = "black", face = "italic", vjust = 0.5, hjust = 0.5)) +
		theme(axis.title.y = element_text(size = 10, color = "black", face = "italic", vjust = 0.5, hjust = 0.5)) +
		theme(axis.text.x = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
		theme(axis.text.y = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
		theme(legend.text = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
	)
}
dev.off()
