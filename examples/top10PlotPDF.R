#!/usr/bin/Rscript --vanilla

# (c) 2015 Charite Universtitaetsmedizin Berlin
# Author:  Leon Kuchenbecker <lkuchenb@inf.fu-berlin.de>
# License: GPLv2

require(ggplot2, quietly = T)
require(reshape2, quietly = T)

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 2) stop("Usage: top10PlotPDF <input file> [<input file> [...]] <output file>")

outputPath <- rev(argv)[1]
argv <- head(argv, -1)

if (file.exists(outputPath)) stop("Output file exists!")

readCountFile <- function(path) {
  d <- read.table(path)
  d <- cbind(do.call(rbind.data.frame, strsplit(as.character(d[,1]), ":")), d[,2])
  colnames(d) <- c("V","cdrAA","J", path)
  d[,4] <- d[,4] / sum(d[,4])
  return(d)  
}

# Read and merge data
data <- readCountFile(argv[1])
for (path in tail(argv,-1)) {
  data <- merge(data, readCountFile(path), by=c("V","cdrAA","J"), all=T)
}
data <- melt(data, id.vars =  1:3, variable.name = "File", value.name = "Frequency")
data$Frequency[is.na(data$Frequency)] <- 0
data <- data[order(-data$Frequency),]


# Identify top 10 clonotypes in the sample that was specified first at the
# command line
top <- data[data$File==argv[1],c("V","cdrAA","J")][1:10,]
top$Clonotype <- paste0("(",top$V,") ",top$cdrAA," (", top$J,")")
top$Clonotype <- factor(top$Clonotype, levels=unique(top$Clonotype))

top <- merge(top, data, all.x=T)
pdf(outputPath)
ggplot(top) + geom_bar(aes(x=Clonotype, y=Frequency, fill=File), stat="identity", position="dodge") + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_brewer(palette="Set1")
print(last_plot())
invisible(dev.off())

cat("Wrote output to '", outputPath, "'.\n", sep = "")
