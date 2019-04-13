setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("TickDynamics.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varible

#setting up each factor so not treated as a continous variable
data$Instar<-factor(data$Instar)
data$year<-factor(data$year)
data$dayOyear<-factor(data$dayOyear)
data$Date_Collected<-factor(data$Date_Collected)

#subset data
dataIN1<-subset(data,Instar==1)
dataIN2<-subset(data,Instar==2)
dataIN3<-subset(data,Instar==3)

library(MASS)
m1 <- glm.nb(Count ~ Host, data = dataIN1)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~ year, data = dataIN1)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~ Host, data = dataIN2)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~ year, data = dataIN2)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~ Host, data = dataIN3)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~ year, data = dataIN3)
summary(m1)
anova(m1)

m1 <- glm.nb(Count ~  year + X2week , data = dataIN3)
summary(m1)
anova(m1)
Anova(m1,type="II")
