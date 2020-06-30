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
library(mgcv)
library(plyr)
library(boot)

#This is based on a previous GLM

#Run the GAM, with smoother applied to s()
gam1<-gam(co3~s(Jday.1), family = binomial, data=data)
summary(gam1)

#Visualization of the GAM
##Quick look
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
##Pulling out the data
plotdata <- visreg(gam1, type = "contrast", plot = F) 
##Pulling out just the smoother
smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))
##Placing the smoothers in a sinlge variable
P5 <- smooths[ which(smooths$Variable== "Jday.1"), ]
##Correction for BIONOMIAL data
P5$smooth2 <- inv.logit(P5$smooth)
P5$lower2 <- inv.logit(P5$lower)
P5$upper2 <- inv.logit(P5$upper)
##plot
plot(P5$smooth)
lines(P5$upper)
plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)
