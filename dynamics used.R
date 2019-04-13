#this is the actual script used

setwd("D:/Rfiles/TickPrev") #how to set my working directory
Ldata<-read.csv("Adults.csv",header=TRUE) #how to import a file, csv file, and has a header
summary(Ldata) #gives means etc

#setting up each factor so not treated as a continous variable
Ldata$Year<-factor(Ldata$Year)

library(MASS)
library(car)
#needed to extract the AIC to determine best model

m1 <- glm.nb(Count ~ Gen_Location+Host+Year, data = Ldata)
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

m1a<-glm.nb(Count ~ Gen_Location+Host+Year, data = Ldata)
m1b<-glm.nb(Count ~ Host+Year+Gen_Location, data = Ldata)
m1c<-glm.nb(Count ~ Year+Gen_Location+Host, data = Ldata)
Anova(m1b,type = 2)
anova(m1a)

#to get post hoc comparisons
library(multcomp)
summary(glht(m2, mcp(Host="Tukey")))

#to build a general additive model
library(mgcv)
gam1<-gam(Count~s(Day)+Year+Host+Gen_Location, family = negbin(theta= c(1,10)), data=Ldata)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)


#trying something else for graphing gam1 which I did use!
library(visreg)
library(plyr)
plotdata <- visreg(gam1, type = "contrast", plot = F)

smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5 <- smooths[ which(smooths$Variable== "Day"), ]
P5$smooth2 <- exp(P5$smooth)
P5$lower2 <- exp(P5$lower)
P5$upper2 <- exp(P5$upper)

plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)

#graphing with lines
library(ggplot2)
p<-ggplot(Ldata,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limit=c(0.15,120))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()
p+scale_shape_manual(values = c(16:19))

#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(ggplot2)
plot1 <- ggplot(aes(y = (Count), x = factor(Host)), data = Ndata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3,aes(color = Host,shape=Year)) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10()+
  theme_bw()
plot1
plot1+scale_shape_manual(values = c(15:19))

##Nyphs data
Ndata<-read.csv("Nymphs.csv",header=TRUE) #how to import a file, csv file, and has a header
summary(Ndata) #gives means etc

#setting up each factor so not treated as a continous variable
Ndata$Year<-factor(Ndata$Year)

library(MASS)
library(car)
#needed to extract the AIC to determine best model

m1 <- glm.nb(Count ~ Gen_Location + Year +Host, data = Ndata)
m2 <- glm.nb(Count ~ Gen_Location, data = Ndata)
m3 <- glm.nb(Count ~ Gen_Location+Host, data = Ndata)
m4 <- glm.nb(Count ~ Host+Year, data = Ndata)
summary(m1)
AIC(m3)
Anova(m4,type="III")
anova(m4)

#to get post hoc comparisons
library(multcomp)
summary(glht(m4, mcp(Year="Tukey")))

#to build a general additive model
library(mgcv)
gam1<-gam(Count~s(Day)+Year+Host, family = negbin(theta= c(1,10)), data=Ndata)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)


#trying something else for graphing gam1 which I did use!
library(visreg)
plotdata <- visreg(gam1, type = "contrast", plot = F)

smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5 <- smooths[ which(smooths$Variable== "Day"), ]
P5$smooth2 <- exp(P5$smooth)
P5$lower2 <- exp(P5$lower)
P5$upper2 <- exp(P5$upper)

plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)

#graphing with lines
ggplot(Ndata,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limits=c(min(P5$smooth2),max(Ndata$Count)))+
  scale_x_continuous(limits=c(150,220),breaks = c(151,181,212))+
  theme_bw()

#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(ggplot2)
plot1 <- ggplot(aes(y = (Count), x = factor(Year)), data = Ndata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10(limits=c(min(P5$smooth2),max(Ndata$Count)))+
  theme_bw()
plot1


##Adult data
Adata<-read.csv("Adults.csv",header=TRUE) #how to import a file, csv file, and has a header
summary(Adata) #gives means etc

#setting up each factor so not treated as a continous variable
Adata$Year<-factor(Adata$Year)

library(MASS)
library(car)
#needed to extract the AIC to determine best model

m1 <- glm.nb(Count ~ Gen_Location + Year +Host, data = Adata)
m2 <- glm.nb(Count ~ Host, data = Adata)
m3 <- glm.nb(Count ~ Year, data = Adata)
m4 <- glm.nb(Count ~ Host+Year, data = Adata)
summary(m1)
  AIC(m2)
Anova(m4,type="III")
anova(m4)

#to get post hoc comparisons
library(multcomp)
summary(glht(m4, mcp(Year="Tukey")))

#to build a general additive model
library(mgcv)
gam1<-gam(Count~s(Day)+Year+Host, family = negbin(theta= c(1,10)), data=Adata)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)


#trying something else for graphing gam1 which I did use!
library(visreg)
plotdata <- visreg(gam1, type = "contrast", plot = F)

smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5 <- smooths[ which(smooths$Variable== "Day"), ]
P5$smooth2 <- exp(P5$smooth)
P5$lower2 <- exp(P5$lower)
P5$upper2 <- exp(P5$upper)

plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)

#graphing with lines
ggplot(Adata,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10()+
  scale_x_continuous(breaks = c(151,181,212))+
  theme_bw()

#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

plot1 <- ggplot(aes(y = (Count), x = factor(Year)), data = Adata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10()+
  theme_bw()
plot1


##Looking at sex ratios

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("AdultsSR2.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varible

#setting up each factor so not treated as a continous variable
data$Year_Collected<-factor(data$Year_Collected)
data$FM<-factor(data$FM)
data$Eng<-factor(data$Eng)

#GLM time
library(MASS)
library(car)
m1 <- with(data, glm(FM ~ Gen_Location+Host+Year_Collected, family = binomial))
m2 <- with(data, glm(FM ~ Gen_Location+Host, family = binomial))
m3 <- with(data, glm(FM ~ Gen_Location+Year_Collected, family = binomial))
m4 <- with(data, glm(FM ~ Host+Year_Collected, family = binomial))
m5 <- with(data, glm(FM ~ Gen_Location, family = binomial))
m6 <- with(data, glm(FM ~ Host, family = binomial))
m7 <- with(data, glm(FM ~ Year_Collected, family = binomial))
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)
logLik(m6)
logLik(m7)
AIC(m1,m2,m3,m4,m5,m6,m7)



library(multcomp)
summary(glht(m3, mcp(Host="Tukey")))


#to build a general additive model
library(mgcv)
gam1<-gam(Eng~s(Day)+Gen_Location+Year_Collected, family = binomial, data=data)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)

library(visreg)
plotdata <- visreg(gam1, type = "contrast", plot = F)

library(plyr)
smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))


library(boot)
P5 <- smooths[ which(smooths$Variable== "Day"), ]
P5$smooth2 <- inv.logit(P5$smooth)
P5$lower2 <- inv.logit(P5$lower)
P5$upper2 <- inv.logit(P5$upper)

plot(P5$smooth)
lines(P5$upper)
plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)

#graphing with lines
library(ggplot2)
#have to reload data so that it is a sex is 1 and0

fig<-ggplot(data,aes(x = Day, y = Eng))+
  geom_line(aes(x, smooth2), P5)+
  ylab("Ratio Female")+
  xlab("Julian Day")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(limit=c(145,230),breaks = c(151,181,212))+
  theme_bw()
fig

fig+stat_summary(fun.data = "mean_cl_boot")

##Looking at engorgment number

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("larvaE.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varible

data<-subset(data,Age=="L")

#setting up each factor so not treated as a continous variable
data$Year_Collected<-factor(data$Year_Collected)
data$Eng<-factor(data$Eng)

#GLM time
library(MASS)
library(car)

m1 <- with(data, glm(Eng ~ Gen_Location+Host+Year_Collected, family = binomial))
m2 <- with(data, glm(Eng ~ Gen_Location+Host, family = binomial))
m3 <- with(data, glm(Eng ~ Gen_Location+Year_Collected, family = binomial))
m4 <- with(data, glm(Eng ~ Host+Year_Collected, family = binomial))
m5 <- with(data, glm(Eng ~ Gen_Location, family = binomial))
m6 <- with(data, glm(Eng ~ Host, family = binomial))
m7 <- with(data, glm(Eng ~ Year_Collected, family = binomial))
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)
logLik(m6)
logLik(m7)
AIC(m1,m2,m3,m4,m5,m6,m7)

m1<-with(data, glm(Eng~Year_Collected, family = binomial))
m2<-with(data, glm(Eng~Host+Year_Collected, family = binomial))
m3<-with(data, glm(Eng~Gen_Location+Year_Collected, family = binomial))
m4<-with(data, glm(Eng~Gen_Location, family = binomial))
AIC(m1,m2,m3,m4)
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
Anova(m1,type=3)

library(multcomp)
summary(glht(m6, mcp(Host="Tukey")))


#to build a general additive model
library(mgcv)
gam1<-gam(Eng~s(Day)+Host+Year_Collected+Gen_Location, family = binomial, data=data)
summary(gam1)
AIC(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)

library(visreg)
plotdata <- visreg(gam1, type = "contrast", plot = F)

library(plyr)
smooths <- ldply(plotdata, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))


library(boot)
P5 <- smooths[ which(smooths$Variable== "Day"), ]
P5$smooth2 <- inv.logit(P5$smooth)
P5$lower2 <- inv.logit(P5$lower)
P5$upper2 <- inv.logit(P5$upper)

plot(P5$smooth)
lines(P5$upper)
plot(P5$x,P5$smooth2)
lines(P5$x,P5$lower2)
lines(P5$x,P5$upper2)

#graphing with lines
library(ggplot2)
#have to reload data so that it is a sex is 1 and0

fig<-ggplot(data,aes(x = Day, y = Eng))+
  geom_line(aes(x, smooth2), P5)+
  ylab("Ratio Female")+
  xlab("Julian Day")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(limit=c(145,230),breaks = c(151,181,212))+
  theme_bw()
fig

fig+stat_summary(fun.data = "mean_cl_boot")

data<-read.csv("adultFE.csv",header=TRUE)
gam1<-gam(Eng~s(Day), family = binomial, data=data)
plotdata <- visreg(gam1)
plotdata$fit$visregFit2<-inv.logit(plotdata$fit$visregFit)


fig<-ggplot(plotdata$fit,aes(x = Day, y = vistregFig))+
  geom_line(aes(plotdata$fit$Day,(plotdata$fit$visregFit2)))+
  ylab("Engorged")+
  xlab("Julian Day")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(limit=c(145,230),breaks = c(150,160,170,180,190,200,212))+
  theme_bw()
fig

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("nymphE.csv",header=TRUE) 
summary(data)
data$Year<-factor(data$Year_Collected)
summary(data)
data<-na.omit(data)

prev<-data[data$Gen_Location=="LFOGO",9] ##must change data name, factor name, and specific factor, 14 refers to pos_neg column
length(prev)
summary(prev)
emptyvar<-numeric(10000)
for (i in 1:10000) {emptyvar[i]<-mean(sample(prev,replace=T),)}
quantile(emptyvar, c(0.025, 0.975))
mean(prev)
mean(emptyvar)
