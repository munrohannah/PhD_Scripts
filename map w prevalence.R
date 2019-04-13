library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency

#drawing large map
map('worldHires','Canada', xlim=c(-61,-50), ylim=c(46,57), fill=TRUE, col="gray90")

#To add points to the file i have lat and long in a file.
setwd("D:/Rfiles/Maps")
data<-read.csv("studysitesprev.csv",header=TRUE)

#Trying to add pie charts
#FAIL