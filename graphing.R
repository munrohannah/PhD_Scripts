#graphing seasonality of just comu in witless bay

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("COMUwit.csv",header=TRUE) #how to import a file, csv file, and has a header
summary(data) #gives means etc

#subsetting data for each instar
dataL<-subset(data,instar==1)
dataN<-subset(data,instar==2)
dataA<-subset(data,instar==3)

library(MASS)
library(car)
library(mgcv)
library(visreg)
library(plyr)
library(ggplot2)
library(boot)

gamL<-gam(Count~s(Day), family = negbin(theta= c(1,10)), data=dataL)
gamN<-gam(Count~s(Day), family = negbin(theta= c(1,10)), data=dataN)
gamA<-gam(Count~s(Day), family = negbin(theta= c(1,10)), data=dataA)

plotdataL <- visreg(gamL)
plotdataN <- visreg(gamN)
plotdataA <- visreg(gamA)

plotdataL$fit$smooth<-exp(plotdataL$fit$visregFit)
plotdataL$fit$up<-exp(plotdataL$fit$visregUpr)
plotdataL$fit$lower<-exp(plotdataL$fit$visregLwr)
plotdataL$res$res<-exp(plotdataL$res$visregRes)
plotdataL$res$Pos<-exp(plotdataL$res$visregPos)

plot(plotdataL$fit$Day,log10(exp(plotdataL$fit$visregFit+20)))
lines(plotdataL$fit$Day,log10(exp(plotdataL$fit$visregUpr+20)))
lines(plotdataL$fit$Day,log10(exp(plotdataL$fit$visregLwr+20)))
lines(dataL$Day,dataL$Max_temp)

P5L$smooth2 <- exp(P5L$smooth+1)
P5N$smooth2 <- exp(P5N$smooth+2)
summary(P5N)
P5N$smooth2 <- P5N$smooth2 + 1
P5A$smooth2 <- (exp(P5A$smooth+1))