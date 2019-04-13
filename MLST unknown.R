#this is the actual script used

setwd("D:/Rfiles/MLST") #how to set my working directory
data<-read.delim("STdist.txt",header=TRUE) #how to import a file, csv file, and has a header
summary(data) #gives means etc

#setting up each factor so not treated as a continous variable
data$year<-factor(data$year)

library(MASS)
library(car)
#glm time

m1 <- glm(old ~ Gen_Location+Host+Year, data = Ldata)
m2 <- glm.nb(Count ~ Gen_Location+Host, data = Adata)
m3 <- glm.nb(Count ~ Gen_Location+Year, data = Ldata)
m4 <- glm.nb(Count ~ Host+Year, data = Ldata)
m5 <- glm.nb(Count ~ Gen_Location, data = Ldata)
m6 <- glm.nb(Count ~ Host, data = Ldata)
m7 <- glm.nb(Count ~ Year, data = Ldata)
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)
logLik(m6)
logLik(m7)
AIC(m1,m2,m3,m4,m5,m6,m7)

Anova(m1,type="2")
anova(m1)
