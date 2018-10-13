#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
fileName = args[1] #strsplit(args[1],"/")
jpegName = args[4] #print(tail(name[[1]], 1))

fileNameLength = nchar(fileName)
N = fileNameLength -10
title = substr(fileName, 4, N)

options(stringAsFactors=F)
d <- read.csv(fileName)
library(limma); library(edgeR)
d <- data.frame(d[2:ncol(d)], row.names = d$X)

#design <- cbind(C=1, TvsC=(c(0,0,0, 1,1,1, 0,0,0, 1,1,1)))
#design <- cbind(C=1, TvsC=(c(0,0,0, 0,0,0, 1,1,1, 1,1,1)))
c <- as.numeric(args[2])
t <- as.numeric(args[3])
design <- cbind(C=1, TvsC=(c(rep(0,c), rep(1,t))))

# normalize
v <- voom(d, design, plot=TRUE)
#v <- cpm(d*10, log=T)

# linear fitting and empirical bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=ncol(design), n=nrow(d))

# get subgroups
positives <- tt[tt$logFC>log(1.5) & tt$adj.P.Val <= 0.05,]
negatives <- tt[tt$logFC < -log(1.5) & tt$adj.P.Val<= 0.05,]

# save data in csv file
csvName = paste("volcano", title, ".csv", sep='')
write.csv(tt, file = csvName)
 
# plot
jpeg(jpegName)
plot(tt$logFC, -log(tt$adj.P.Val), cex=.5, xlim=c(-3,3), pch=19, col=rgb(red=.1, green=.1, blue=.1, alpha=0.5),
     ylab = 'log (adj.P-value)', xlab = 'log (Fold change)', main=title)
abline(h=-log(0.05), col='red', lty=2, lwd=2)
abline(v = 0.4, col='red', lty=2, lwd=2)
abline(v = -0.4, col='red', lty=2, lwd=2)

points(positives$logFC, -log(positives$adj.P.Val), col=rgb(red=0, green=1, blue=0, alpha = 0.5))
points(negatives$logFC, -log(negatives$adj.P.Val), col=rgb(red=1, green=0, blue=0, alpha = 0.5))
text(2, 5, dim(positives)[1], col='green')
text(-2, 5, dim(negatives)[1], col='red')

dev.off()
#savePlot(filename = 'Med2_15_GC_limma.jpg')
