#Be sure to have as working directory the directory with the files
setwd("D:/Rfiles/tutorials") #how to set my working directory

#Let's get the filenames
files <- list.files(pattern = "*_Res.csv", include.dirs = TRUE, recursive = TRUE, full.names = TRUE)
files

#Read those files - IMPORTANT: Use lapply instead of sapply.
sets <- lapply(files, read.table, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
names(sets) <- c("S1", "S2", "S3")
str(sets)

#Check if data was read correctly
lapply(sets, dim)
lapply(sets, head)

#Select only those genes with an absolute logFC > 2
lstVenn <- lapply(sets, function(s){
  s[abs(s[,2])> 2, 1]
})

#Check the list just created
lapply(lstVenn, length)
lapply(lstVenn, head)

#Create a vector containing all genes
universe <- union(union(sets$S1[,1], sets$S2[,1]), sets$S3[,1])
length(universe)

head(universe)

#Load gplots
library(gplots)

#Create the Venn diagram (works for up to 5 sets)
groups <- venn(data = lstVenn, universe = universe)
groups

#Limit diagram to some genes
venn(data = lstVenn, universe = universe[1:100])

#Write diagram to a pdf file
pdf("VennDiagram.pdf")
venn(data = lstVenn, universe = universe)
dev.off()

############### Bubble Chart

data <- read.table("dataBubbleChart.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

str(data)
head(data)

#Scale based on size
radius <- sqrt( data$Size / pi )

#Create the bubble chart
#pdf("BubbleChart.pdf")

symbols(data$Quality, data$Preservation, circles=radius, inches=0.35, fg="white", bg=data$Module, xlab="Quality", ylab="Preservation")
text(data$Quality, data$Preservation,data$Module, cex=0.6)
abline(h = 2, col = "gray", lty = "dashed")
abline(h = 10, col = "gray", lty = "dashed")

abline(v = 2, col = "yellow", lty = "dashed")

#dev.off()



############### Balloon Plot
data2 <- read.table("dataBalloonPlot.csv", sep = "\t", header = TRUE)
str(data2)
head(data2)


colors <- vector(mode = "character", length = length(data2$cpue))
colors <- ifelse(data2$cpue > 10, "red", 
                 ifelse(data2$cpue > 5, "orange", 
                        ifelse(data2$cpue > 2, "yellow", "skyblue")))

head(colors)
library(gplots)
pdf("BalloonPlot.pdf")
balloonplot(data2$managementzone, data2$species, data2$cpue, xlab = "Zone", ylab = "Species", main = "CPUE per Zone and Species", text.size = 0.5, label = TRUE, label.digits = 1, cum.margins = FALSE, label.lines = FALSE, colsrt = 40, dotcolor = colors)
dev.off()