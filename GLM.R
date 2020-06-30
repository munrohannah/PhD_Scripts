#Basic opening of files
setwd("D:/Rfiles/ticks") #working directory
data<-read.csv("genetics.csv",header=TRUE) #how to import a file, csv file, and has a header
##Quick look at the data to see that it imported correctly
summary(data)
head(data)
names(data)
##Changing between numeric and factor
data$Year<-factor(data$Year)

#Loading needed packages
library(MASS)
library(car)

#Writing models for AIC
##A series of Binomial models for AIC analysis
m1<-glm(co3~Island+Year+Host, family=binomial,  data= data )
m2<-glm(co3~Year+Host, family=binomial,  data= data)
m3<-glm(co3~Host+Island, family=binomial,  data= data)
m4<-glm(co3~Island+Year, family=binomial,  data= data)
m5<-glm(co3~Year, family=binomial,  data= data)
m6<-glm(co3~Host, family=binomial,  data= data)
m7<-glm(co3~Island, family=binomial,  data= data)
##Extracting the AIC values
AIC(m1,m2,m3,m4,m5,m6,m7)
##
