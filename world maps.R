#This is an attempt to map where we have sequecine data
rm(list=ls())

library(maps) # Provides functions that let us plot the maps 
library(mapdata) # Contains the hi-resolution points that mark out the countries.

map('worldHires') #along international dateline
map('world2Hires') #along prime merdian

map.scale(160,0,relwidth = 0.15, metric = TRUE, ratio = TRUE) 
map.scale(160,-40,relwidth = 0.15, metric = TRUE, ratio = TRUE)

map('worldHires','Italy')
map('worldHires','ch')

map('worldHires', c('UK', 
                    'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'))
map('worldHires', c('UK', 'Ireland', 'Isle of Man','Isle of Wight'), xlim=c(-11,3), ylim=c(49,60.9))  
points(-1.615672,54.977768,col=2,pch=18)


setwd("D:/Rfiles/Maps")
data<-read.csv("Location.csv",header=TRUE)
data
points(data$Long, data$Lat, pch=19, col="red", cex=0.5) 
