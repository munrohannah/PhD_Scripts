setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("PosIndex.csv",header=T) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

Index1GLM<-with(data, glm(Pos_Neg~Age+Index1, family = binomial))
Index1GLM
anova(Index1GLM)
summary(Index1GLM)
plot(Index1GLM)