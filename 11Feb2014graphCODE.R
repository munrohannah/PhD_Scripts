setwd("D:/Rfiles/qPCR") #how to set my working directory
data<-read.csv("11Feb2014b.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
library(ggplot2) #http://www.cookbook-r.com/ is the best resource for plots
g<-ggplot(data,aes(x=Cycle,y=Rn))+geom_point(shape=1)+facet_grid(Row ~ Col)
g 