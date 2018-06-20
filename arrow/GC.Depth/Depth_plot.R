library(ggplot2)
args <- commandArgs()
inputFile <- args[6]
y <- seq(0,0.3, 0.02)
yScalePercent <- y*100
geneExp <- read.table(inputFile)

p <- ggplot(geneExp, aes(geneExp[,1]))
p <- p + geom_histogram(fill="#60BFB9",colour = "#FFFFFF", alpha=1, binwidth=1,position="identity",aes(y=..density..))
#p <- p + geom_vline(aes(xintercept=mean(geneExp[,1], na.rm=T)), colour="#EC5A6C", linetype="dashed", size=1)
p <- p + theme_bw()
p <- p + scale_x_continuous(limits=c(0,150), breaks=seq(0,150,10))
p <- p + scale_y_continuous(breaks = y, labels = yScalePercent ) 
p <- p + theme(panel.grid.major=element_line(size=0.25,linetype =3,color="#9FA0A0"))
p <- p + theme(panel.grid.minor=element_line(size=0.25,linetype =3,color="#9FA0A0"))
p <- p + theme(panel.border=element_rect(fill='transparent', color='black',size=1))
#p <- p + theme(panel.border=element_line(size=1))
p <- p + theme(plot.title = element_text(size=25))
p <- p + theme(axis.text = element_text(size = 18))
p <- p + theme(axis.title = element_text(size = 20))
p <- p + xlab("Sequencing depth (X)")
p <- p + ylab("Percent of bases (%)")
ggsave(paste(inputFile, ".png", sep=""), width=12, height=8)

