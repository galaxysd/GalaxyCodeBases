library(ggplot2)
library(splines)
args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]

pdf(output)
file.tab <- read.table(file=input,sep ="\t",header = F)
file <- file.tab[,1]

for (i in 1:length(file)){
	infile = as.character(file[i])
	head <- unlist(strsplit(infile, split = "[/.]"))
	input.tab <- try(read.table(file=infile,sep ="\t",header = F))
	if(class(input.tab)=='try-error'){
		next
	}
	print(
		ggplot(data=input.tab,aes(x=input.tab[,1],y=input.tab[,2])) +
#		geom_histogram(position="dodge", bins = 200,na.rm = TRUE, color="darkblue", fill="lightblue") +
		geom_bar(stat="identity",position="dodge",na.rm = TRUE, color="darkblue", fill="lightblue") +
		labs(title=head[length(head)-1],x="Length",y="Number of Fragment",size=20) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#		scale_y_continuous(breaks=seq(0,10000,100)) +
		scale_x_continuous(breaks=seq(0,1000,100)) +
		theme(axis.title.x = element_text(size = 10, color = "black", face = "italic", vjust = 0.5, hjust = 0.5)) +
		theme(axis.title.y = element_text(size = 10, color = "black", face = "italic", vjust = 0.5, hjust = 0.5)) +
		theme(axis.text.x = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
		theme(axis.text.y = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
		theme(legend.text = element_text(size = 8, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
	)
}
dev.off()
