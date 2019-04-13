setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("Nov3data.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

#setting up each factor so not treated as a continous variable
data$A.E<-factor(data$A.E)
data$Month<-factor(data$Month)
data$Year<-factor(data$Year)
data$Extraction_Date<-factor(data$Extraction_Date)

#remove the AE=1 for analysis because there are only two observations
dataAE2.7<-subset(data,A.E!=1)
head(dataAE2.7)
str(dataAE2.7)

#gam to see the influence of date
library(mgcv)
gam1<-gam(pos_neg~s(Jday)+Year+A.E+Host+Gen_Location, family = binomial, data=dataAE2.7)
summary(gam1)
plot(gam1)

#placing summary in cdata
#summarizing the data based on differnet factors.
library(plyr)
library(MASS)
cdata <- ddply(dataAE2.7, c("A.E") , summarise, #just need to change the factor
               N    = length(pos_neg),
               mean = mean(pos_neg),
               count = sum(pos_neg)
               )

#calaculating the 95CI using the wald formula. not perfect if you have skewed data
library(prevalence)
CI<-propCI(cdata$count,cdata$N,method="wald")##calculate the CI for each in cdata
cdata$upper<-CI$upper #placing the upper 95CI in cdata
cdata$lower<-CI$lower

AEdata<-cdata #renaming

#do not use this, does not work on skewed data, should be similar to the wald calaculations above
cdata$ME<-qnorm(0.975)*sqrt(cdata$mean*(1-cdata$mean)/cdata$N) #margine of error
cdata$lower95<-cdata$mean-cdata$ME #lower limit
cdata$upper95<-cdata$mean+cdata$ME #upper limit
cdata


#pulling out the pos_neg data based on differnet factors.
AE2<-data[data$A.E==2,14] ##must change data name, factor name, and specific factor, 14 refers to pos_neg column

AE2<-data[data$A.E==2,14]
AE3<-data[data$A.E==3,14]
AE4<-data[data$A.E==4,14]
AE5<-data[data$A.E==5,14]
AE6<-data[data$A.E==6,14]
AE7<-data[data$A.E==7,14]

HostATPU<-data[data$Host== "ATPU",14]
HostBLKI<-data[data$Host== "BLKI" ,14]
HostCOMU<-data[data$Host== "COMU" ,14]
HostRAZO<-data[data$Host== "RAZO" ,14]
HostTBMU<-data[data$Host== "TBMU" ,14]
HostUnknown<-data[data$Host== "Unknown" ,14]

Month6<-data[data$Month== "6",14]
Month7<-data[data$Month== "7" ,14]
Month8<-data[data$Month== "8" ,14]

Year2011<-data[data$Year== "2011",14]
Year2012<-data[data$Year== "2012" ,14]
Year2013<-data[data$Year== "2013" ,14]
Year2014<-data[data$Year== "2014" ,14]

Gen_LocationGANN<-data[data$Gen_Location== "GANN",14]
Gen_LocationGREAT<-data[data$Gen_Location== "GREAT" ,14]
Gen_LocationGULL<-data[data$Gen_Location== "GULL" ,14]
Gen_LocationLFOGO<-data[data$Gen_Location== "LFOGO" ,14]

ALL<-data[,14]

#alex's suggestion. does not work perfectly as it just pulls out 0s and 1
fake.data <- sample(AE2, 100000, replace=TRUE, prob=NULL)
quantile(fake.data, c(0.025, 0.975))


#randomization calcluation of 95CI
emptyvar<-numeric(1000)                     #build an empty variable for the prevalence to be placed in
for (i in 1:1000) {                         #starting a loop
  emptyvar[i]<-mean(sample(AE2,replace=T),) #place mean of randomly sampled data, must replace vector
}                                           #close the loop
AE2CI<-quantile(emptyvar, c(0.025, 0.975))   #give the 95%CI

#AE
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE2,replace=T),)}
AE2CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE3,replace=T),)}
AE3CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE4,replace=T),)}
AE4CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE5,replace=T),)}
AE5CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE6,replace=T),)}
AE6CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(AE7,replace=T),)}
AE7CI<-quantile(emptyvar, c(0.025, 0.975))

#Host
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostATPU,replace=T),)}
HostATPUCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostBLKI,replace=T),)}
HostBLKICI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostCOMU,replace=T),)}
HostCOMUCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostRAZO,replace=T),)}
HostRAZOCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostTBMU,replace=T),)}
HostTBMUCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(HostUnknown,replace=T),)}
HostUknownCI<-quantile(emptyvar, c(0.025, 0.975))

#Month
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Month6,replace=T),)}
Month6CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Month7,replace=T),)}
Month7CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Month8,replace=T),)}
Month8CI<-quantile(emptyvar, c(0.025, 0.975))

#Year
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Year2011,replace=T),)}
Year2011CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Year2012,replace=T),)}
Year2012CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Year2013,replace=T),)}
Year2013CI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Year2014,replace=T),)}
Year2014CI<-quantile(emptyvar, c(0.025, 0.975))

#Gen_Location
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Gen_LocationGANN,replace=T),)}
Gen_LocationGANNCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Gen_LocationGREAT,replace=T),)}
Gen_LocationGREATCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Gen_LocationGULL,replace=T),)}
Gen_LocationGULLCI<-quantile(emptyvar, c(0.025, 0.975))
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(Gen_LocationLFOGO,replace=T),)}
Gen_LocationLFOGOCI<-quantile(emptyvar, c(0.025, 0.975))

emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(ALL,replace=T),)}
ALLCI<-quantile(emptyvar, c(0.025, 0.975))
ALLCI


#appending CI to the end of the cdata dataframe
AECI<-data.frame(t(data.frame(AE2CI,AE3CI,AE4CI,AE5CI,AE6CI,AE7CI)))
cdata$CI2.5<-AECI$X2.5.
cdata$CI97.5<-AECI$X97.5.

HostCI<-data.frame(t(data.frame(HostATPUCI,HostBLKICI,HostCOMUCI,HostRAZOCI,HostTBMUCI,HostUknownCI)))
cdata$CI2.5<-HostCI$X2.5.
cdata$CI97.5<-HostCI$X97.5.

MonthCI<-data.frame(t(data.frame(Month6CI,Month7CI,Month8CI)))
cdata$CI2.5<-MonthCI$X2.5.
cdata$CI97.5<-MonthCI$X97.5.

YearCI<-data.frame(t(data.frame(Year2011CI,Year2012CI,Year2013CI,Year2014CI)))
cdata$CI2.5<-YearCI$X2.5.
cdata$CI97.5<-YearCI$X97.5.

Gen_LocationCI<-data.frame(t(data.frame(Gen_LocationGANNCI,Gen_LocationGREATCI,Gen_LocationGULLCI,Gen_LocationLFOGOCI)))
cdata$CI2.5<-Gen_LocationCI$X2.5.
cdata$CI97.5<-Gen_LocationCI$X97.5.

Gen_Locationdata<-cdata #renaming

#graphing time
library(ggplot2)

p <- ggplot(Monthdata, aes(Month, mean))#for single factors
p + geom_pointrange(aes(ymin = CI2.5, ymax = CI97.5))+ ylim(0,0.4) + theme_minimal()

ggsave("MonthB.svg",width=5,height=5)

p <- ggplot(Yeardata, aes(Year, mean))#for single factors
p + geom_crossbar(aes(ymin = CI2.5, ymax = CI97.5),width = 0.4)+ expand_limits(y=0,y=0.4) + theme_minimal()

p <- ggplot(cdata, aes(Host, mean, colour = A.E))#for multiple factors
p + geom_crossbar(aes(ymin = lower95, ymax = upper95),width = 0.4, position = "dodge")

write.csv(cdata,file = "cdata.csv")
