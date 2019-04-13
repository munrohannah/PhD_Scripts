setwd("D:/Rfiles/AI") #how to set my working directory
data<-read.csv("RawGull09.11.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
model1<-glm( Pos.Neg ~ Species+Year+Season, data=data, family = binomial)
model1
summary(model1)
anova(model1)