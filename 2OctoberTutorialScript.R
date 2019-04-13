setwd("D:/Rfiles/tutorials") #how to set my working directory
flowers<-read.csv("flowers.csv",header=TRUE) #how to import a file, csv file, and has a header
head(flowers) # gives the first few lines of data
str(flowers) #check the data type and varibles
test1<-lm(flowers$Sepal.Length~flowers$Species)
test2<-lm(Sepal.Length~Species,data=flowers)
plot(test2)

plot(Sepal.Length~Species, data=flowers)