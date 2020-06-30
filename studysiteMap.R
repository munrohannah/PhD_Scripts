#from http://www.molecularecologist.com/2012/09/making-maps-with-r/
library(maps) #Maping package
library(mapdata) #Place with the canada map
library(maptools)
library (scales)

#I am going to draw a map here. The X is your long and y is your lat in decimals. the col is the color.
map('worldHires','Canada', xlim=c(-52.9,-52.7), ylim=c(47.1,47.3), fill=TRUE, col="gray90")

#To add points to the file i have lat and long in a file.
setwd("D:/Rfiles/Maps")
data<-read.csv("studysitesWB.csv",header=TRUE)
data
points(data$Long, data$Lat, pch=19, col="red", cex=0.5)  #plot my sample sites

#Not to at the lables. the offset is how far the title is from the point.
library(calibrate)
textxy(data$Long, data$Lat, labs=data$Title, offset=-0.5)


#For the smaller map
library(raster)
# Finding ISO3 code for Canada
getData('ISO3')  # Canada's code is "CAN"
# Reading in data at different levels
can0<-getData('GADM', country="CAN", level=0)
can1<-getData('GADM', country="CAN", level=1)
can2<-getData('GADM', country="CAN", level=2)
class(can0)
str(can0)
class(can1)
str(can1)
plot(can0,xlim=c(-52.8,-52.7), ylim=c(47.15,47.3))

#To add points to the file i have lat and long in a file.
setwd("D:/Rfiles/Maps")
data<-read.csv("studysitesWB.csv",header=TRUE)
data
points(data$Long, data$Lat, pch=19, col="red", cex=0.5)  #plot my sample sites

#Not to at the lables. the offset is how far the title is from the point.
library(calibrate)
textxy(data$Long, data$Lat, labs=data$Title, offset=-0.5)
