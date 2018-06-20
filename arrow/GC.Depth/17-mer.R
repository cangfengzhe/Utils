## read kmer.freq file
args <- commandArgs(trailingOnly=TRUE)
file = args[1]		#DCXC.freq.stat
cw=getwd()  #get current work directory
#kmerdata<-read.table(file = "17mer.freq",list(A="Depth",B="Frequence",C="Frequency",D="Depth*Frequency"), 
#    skip = 0, nlines = 255,na.strings = "NA",blank.lines.skip = TRUE, multi.line = TRUE)
data <- read.table(file, skip=8)#, nlines=255, na.strings="NA", blank.lines.skip=TRUE, multi.line = TRUE)
#data <- as.matrix(data)
data[,3]
kmer <- data[,3] * 100
kmer
#data <- as.matrix(data)

## drawing
xrange  <- 200
yrange  <- 5

pdf("17-mer.pdf", width=9, height=6)
par(font=2, font.lab=2, font.axis=1.5, lwd=1.5)

#plot(kmerdata$A,kmerdata$C,type="l",col="blue",xlab="",ylab="",cex=0.25,pch=18,xlim=c(0,xrange),ylim=c(0,yrange),
#     main="17-mer Depth Distribution Curve",lwd=4,xaxt="n",yaxt="n",cex.main=1.3)
plot(#c(1:length(kmer)), 
     kmer, 
     type="l", 
     col="blue", 
     xlab="",   ylab="",
     cex=0.25,
     pch=18,
     xlim=c(0,xrange), ylim=c(0,yrange),
     main="17-mer Depth Distribution Curve",
     lwd=2,
     xaxt="n",  yaxt="n",
     cex.main=1.3
     )
xlabel=seq(0,xrange,by=5)
axis(1,at=xlabel,labels=xlabel,las=1,cex.axis=1)
axis(2,las=1,seq(0,yrange,by=0.5),cex.axis=1)
mtext("Frequence(%)",line=2.4,side=2,cex=1.1)
mtext("Depth(X)",line=2.4,side=1,cex=1.1)
box()
dev.off();
