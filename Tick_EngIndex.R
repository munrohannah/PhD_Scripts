Tick_EngIndex.csv

setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("Tick_EngIndex.csv",header=TRUE)#how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

maleAdult <- data[ which(data$Age=='AM'), ]

library(ggplot2)
p <- ggplot(maleAdult, aes(factor(Pos_Neg), Index1 ))
p+geom_boxplot(aes(fill = factor(Age)))

p + geom_boxplot()

a<-with(data, glm(Pos_Neg~Gen_Location,binomial))
anova(a)
AIC(a)
summary(a)
