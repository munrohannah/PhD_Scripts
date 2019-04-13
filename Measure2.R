setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("measure2_rAnaysis.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
library(ggplot2) #http://www.cookbook-r.com/ is the best resource for plots
#Plot below graphs measurement ID on the x, value on the y. 
#There is differnet shaped markers for Engorged/Not
#The differnet colors are for each person who measured.
g<-ggplot(data,aes(x=Measurement_ID,y=Measurement_value))+geom_point(aes(colour = Measured_by))+aes(shape = factor(Engorged))+ facet_grid(Measurement~Age,scales = "free")
g 
