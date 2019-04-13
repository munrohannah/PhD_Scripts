setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("Nov3data.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

Index1GLM<-with(data, glm(pos_neg~A.E, family = binomial))
Index1GLM
anova(Index1GLM)
summary(Index1GLM)
aov(Index1GLM)
plot(Index1GLM)

Index1GLM<-with(data, glm(pos_neg~Gen_Location+Month+Year+A.E+Host, family = binomial))
Index1GLM<-with(data, glm(pos_neg~Host-1, family = binomial))
car::Anova(Index1GLM, test="LR", type="III")

confint(Index1GLM)

data$A.E<-factor(data$A.E)

Gen_Locationdata

dataAE2.7<-subset(data, A.E >1 )
dataAE2.7$A.E<-factor(dataAE2.7$A.E)

Index1GLM<-with(dataAE2.7, glm(pos_neg~A.E, family = binomial))
summary(Index1GLM)
car::Anova(Index1GLM, test="LR", type="III")

Index1GLM<-with(data, glm(pos_neg~Host, family = binomial))
summary(Index1GLM)
car::Anova(Index1GLM, test="LR", type="III")

Yearsub <- subset(data, Year == c("2013","2014"))#pull out subset of data
AEsub<-subset(data, A.E == c(2,5))

Index1GLM<-with(AEsub, glm(pos_neg~A.E, family = binomial))
summary(Index1GLM)
car::Anova(Index1GLM, test="LR", type="III")
Anova(Index1GLM)
