setwd("D:/Rfiles/qPCR") #how to set my working directory
data<-read.csv("17JanRaw.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
library(ggplot2) #http://www.cookbook-r.com/ is the best resource for plots
g<-ggplot(data,aes(x=Cycle,y=Rn)) #setting up the x and y axies
         +geom_point(shape=1) #setting the point style, must be in there
         +facet_grid(Row ~ Col) #so that it is in a grid
g #to visualize the plot

g<-ggplot(data,aes(x=Cycle,y=Rn))+geom_point(shape=1)+facet_grid(Row ~ Col)
g
