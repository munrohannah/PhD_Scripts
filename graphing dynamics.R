#graphing

#this is the actual script used

setwd("D:/Rfiles/TickPrev") #how to set my working directory
Ldata<-read.csv("Larva.csv",header=TRUE) #how to import a file, csv file, and has a header
Ndata<-read.csv("Nymphs.csv",header=TRUE) #how to import a file, csv file, and has a header
Adata<-read.csv("Adults.csv",header=TRUE) #how to import a file, csv file, and has a header
Alldata<-read.csv("allDyn.csv",header=TRUE) #how to import a file, csv file, and has a header


Ldata$Year<-factor(Ldata$Year)
Ndata$Year<-factor(Ndata$Year)
Adata$Year<-factor(Adata$Year)
Alldata$instar<-factor(Alldata$instar)


library(MASS)
library(car)
library(mgcv)
library(visreg)
library(plyr)
library(ggplot2)
library(boot)

gamL<-gam(Count~s(Day)+Year+Host, family = nb(), data=dataL)
gamN<-gam(Count~s(Day)+Year+Host+PrecipAanom, family = nb(), data=dataN)
gamA<-gam(Count~s(Day)+Host+Loc+PrecipManom, family = nb(), data=dataA)

summary(gamA)
plot(gamA)
plotdataL <- visreg(gamL, type = "contrast", plot = F)
plotdataN <- visreg(gamN, type = "contrast", plot = F)
plotdataA <- visreg(gamA, type = "contrast", plot = F)

smoothsL <- ldply(plotdataL, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

smoothsN <- ldply(plotdataN, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

smoothsA <- ldply(plotdataA, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5L <- smoothsL[ which(smoothsL$Variable== "Day"), ]
P5N <- smoothsN[ which(smoothsN$Variable== "Day"), ]
P5A <- smoothsA[ which(smoothsA$Variable== "Day"), ]

P5L$smooth2 <- exp(P5L$smooth+1)
P5N$smooth2 <- exp(P5N$smooth+2)
summary(P5N)
P5N$smooth2 <- P5N$smooth2 + 1
P5A$smooth2 <- (exp(P5A$smooth+1))

fig<-ggplot(Alldata,aes(x = Day, y = Count))+
  geom_point(aes(shape = instar))+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limit=c(0.4,max(Alldata$Count)),breaks = c(1,10,100))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()+
  scale_shape_manual(values = c(1:3)) 
fig+ geom_line(data=P5A,aes(x=x,y=smooth2),lty=5,lwd=1)+
  geom_line(data=P5N,aes(x=x,y=smooth2),lty=3,lwd=1)+
  geom_line(data=P5L,aes(x=x,y=smooth2),lwd=1)


#Go in and change the L N A
p<-ggplot(Ndata,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5N)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limit=c(min(P5N$smooth2),max(Ndata$Count)),breaks = c(1,10,100))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()
p


#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

plot1 <- ggplot(aes(y = (Count), x = factor(Year)), data = Adata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3,aes(color = Host,shape=Year)) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10(limit=c(min(P5A$smooth2),max(Adata$Count)))+
  theme_bw()
plot1

plot1 <- ggplot(aes(y = (Count), x = factor(Host)), data = Ldata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ylab("Count")+
  xlab("Host")+
  scale_y_log10()+
  theme_bw()
plot1

##Looking at engorged level
LEdata<-read.csv("larvaE.csv",header=TRUE) #how to import a file, csv file, and has a header
NEdata<-read.csv("nymphE.csv",header=TRUE) #how to import a file, csv file, and has a header
AEdata<-read.csv("adultFE.csv",header=TRUE) #how to import a file, csv file, and has a header

LEdata$Eng<-factor(LEdata$Eng)

fig2<-ggplot(LEdata,aes(x=Host, y=Eng))+
  stat_summary(fun.data = 'mean_cl_boot')+
  ylab("Percent engorged")+
  xlab("Year_Collected")+
  scale_y_continuous(limit=c(0,1),breaks = c(0,0.5,1))+
  theme_bw()
fig2

gamL<-gam(Eng~s(Day)+Host, family = binomial, data=LEdata)
gamN<-gam(Eng~s(Day)+Year_Collected+Host, family = binomial, data=NEdata)
gamA<-gam(Eng~s(Day)+Year_Collected+Host+Gen_Location, family = binomial, data=AEdata)

plotdataL <- visreg(gamL, type = "contrast", plot = F)
plotdataN <- visreg(gamN, type = "contrast", plot = F)
plotdataA <- visreg(gamA, type = "contrast", plot = F)

smoothsL <- ldply(plotdataL, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

smoothsN <- ldply(plotdataN, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

smoothsA <- ldply(plotdataA, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5L <- smoothsL[ which(smoothsL$Variable== "Day"), ]
P5N <- smoothsN[ which(smoothsN$Variable== "Day"), ]
P5A <- smoothsA[ which(smoothsA$Variable== "Day"), ]

P5L$smooth2 <- inv.logit(P5L$smooth)
P5N$smooth2 <- inv.logit(P5N$smooth)
P5A$smooth2 <- inv.logit(P5A$smooth)

AllE<-read.csv("allE.csv",header=TRUE)
AllE$Age<-as.factor(AllE$Age)
fig<-ggplot(AllE,aes(x = Day, y = Eng))+
  #geom_point(aes(shape = Age))+
  ylab("Count")+
  xlab("Julian Day")+
  #scale_y_log10(limit=c(0.4,max(AllE$Count)),breaks = c(1,10,100))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()+
  scale_shape_manual(values = c(1:3))+
  stat_summary(fun.data = "mean_cl_boot",aes(shape=Age))
fig+ geom_line(data=P5A,aes(x=x,y=smooth2),lty=5,lwd=0.75,col=1)+
  geom_line(data=P5N,aes(x=x,y=smooth2),lty=3,lwd=0.75,col=10)+
  geom_line(data=P5L,aes(x=x,y=smooth2),lwd=0.75,col=20)
fig


##Looking at sex ratios

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("AdultsSR.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varible

#setting up each factor so not treated as a continous variable
data$Year_Collected<-factor(data$Year_Collected)
data$FM<-factor(data$FM)

#GLM time
library(MASS)
library(car)

m1<-with(data, glm(FM~Gen_Location+Host+Year_Collected, family = binomial))
m2<-with(data, glm(FM~Host, family = binomial))
m3<-with(data, glm(FM~Gen_Location+Host, family = binomial))
m4<-with(data, glm(FM~Gen_Location*Host, family = binomial))
AIC(m1,m2,m3,m4)
Anova(m4,type=3)

library(multcomp)
summary(glht(m3, mcp(Host="Tukey")))


#to build a general additive model
library(mgcv)
gam1<-gam(FM~s(Day)+Host, family = binomial, data=data)
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
data<-read.csv("AdultsSR.csv",header=TRUE)

fig<-ggplot(data,aes(x = Day, y = FM))+
  geom_line(aes(x, smooth2), P5)+
  ylab("Ratio Female")+
  xlab("Julian Day")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(limit=c(145,230),breaks = c(151,181,212))+
  theme_bw()
fig

fig+stat_summary(fun.data = "mean_cl_boot")
fig+stat_smooth()

fig2<-ggplot(data,aes(x=Year_Collected, y=Eng))+
  stat_summary(fun.data = 'mean_cl_boot')+
  ylab("Ratio Female")+
  xlab("Year")+
  scale_y_continuous(limit=c(0,1),breaks = c(0,0.5,1))+
  theme_bw()
fig2

fig+stat_summary(fun.data = "mean_cl_boot",geom = "boxplot")
fig+stat_smooth()
