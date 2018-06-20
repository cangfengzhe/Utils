args <- commandArgs(trailingOnly=TRUE)
file = args[1]
cw=getwd()

df <- read.table(file, header = TRUE)
xrange  <- 80
yrange  <- 150

dcols <- densCols(df, colramp=colorRampPalette(c("black", "white")), nbin = 1000)
df$dens <- col2rgb(dcols)[1,] + 1L
cols <- colorRampPalette(c("RoyalBlue", "orange", "red"), space = "Lab")(256)
df$col <- cols[df$dens]
pdf("GC.Depth.pdf", width=9, height=9)
plot(Avgdepth ~ GCpercent, data=df[order(df$dens),], col=col, ylab="Average depth (X)", xlab="GC content (%)", cex.lab = 1.4, cex.axis = 1.3, pch = 20, ylim = c(0,yrange), xlim = c(10,xrange))
dev.off()
