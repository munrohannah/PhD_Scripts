setwd("D:/Rfiles/16s") #how to set my working directory
data<-read.csv("alpha_joost.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

library(car)
Alpha<-with(data, glm(genus~Loc))
Alpha
Anova(Alpha, type="III")
summary(Alpha)
aov(Alpha)
plot(Alpha)


p <- ggplot(data, aes(Loc, genus))
p + geom_boxplot()
