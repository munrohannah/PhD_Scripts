#this is the actual script used

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("TickDynamicsWanom2.csv",header=TRUE) #how to import a file, csv file, and has a header
summary(data) #gives means etc

#setting up each factor so not treated as a continous variable
data$Year<-factor(data$Year)

#subsetting data for each instar
dataL<-subset(data,instar==1)
dataN<-subset(data,instar==2)
dataA<-subset(data,instar==3)

summary(dataN)



library(MASS)
library(car)
#needed to extract the AIC to determine best model

m1 <- glm.nb(Count ~ Host+PrecipManom, data = dataL)
m2 <- glm.nb(Count ~ Host+PrecipJnanom, data = dataL)
m3 <- glm.nb(Count ~ Host+PrecipJlanom, data = dataL)
m4 <- glm.nb(Count ~ Host+PrecipAanom, data = dataL)
m5 <- glm.nb(Count ~ Year+Host+PrecipManom, data = dataL)
m6 <- glm.nb(Count ~ Year+Host, data = dataL)
m7 <- glm.nb(Count ~ Host+PrecipManom, data = dataL)
AIC(m1,m2,m3,m4,m5,m6,m7)
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)

m5 <- glm.nb(Count ~ Loc+Host+PrecipManom, data = dataA)
m5 <- glm.nb(Count ~ Loc+Host, data = dataA)
m5 <- glm.nb(Count ~ Loc+PrecipManom, data = dataA)
AIC(m5)
Anova(m5,type="3")
anova(m5)

summary(m5)
m5

#to get post hoc comparisons
library(multcomp)
summary(dataA)
summary(glht(m5,mcp(Loc ="Tukey")))
m5
glht(m5,mcp(Loc ="Tukey"))

#to build a general additive model

library(mgcv)
gam1<-gam(Count~s(Day)+Host+Loc+PrecipAanom, family = negbin(theta= c(1,10)), data=dataN)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=3,all.terms=TRUE)
AIC(gam1,gam2)

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
p<-ggplot(dataN,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limit=c(0.15,120))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()
p+scale_shape_manual(values = c(16:19))




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
p<-ggplot(dataN,aes(x = Day, y = Count))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  geom_line(aes(x, smooth2), P5)+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(limit=c(0.15,120))+
  scale_x_continuous(limit=c(145,220),breaks = c(151,181,212))+
  theme_bw()
p
p+scale_shape_manual(values = c(16:19))
stat_summary(data=dataA,fun.data = "mean_cl_boot")

#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

counts<-data[dataL$Year==2011,"Count"] ##must change data name, factor name, and specific factor, 14 refers to pos_neg column
min.mean.sd.max(counts)


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

ggplot(dataA,aes(x = PrecipManom, y = Count))+
  geom_point(aes(color = Host, shape=Loc),size= 3)+
  ylab("Count")+
  xlab("May Precip Anomaly")+
  scale_y_log10()+
  theme_bw()

#box plots for the catagorical
#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

plot1 <- ggplot(aes(y = (Count), x = factor(Year)), data = Ldata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10()+
  theme_bw()
plot1


setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("AdultsSR2W.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varible

#setting up each factor so not treated as a continous variable
data$Year_Collected<-factor(data$Year_Collected)
data$FM<-factor(data$FM)
data$Eng<-factor(data$Eng)

#GLM time
library(MASS)
library(car)
m1 <- with(data, glm(FM ~ Gen_Location+Host+PrecipManom, family = binomial))
m2 <- with(data, glm(FM ~ Gen_Location+Host+PrecipJnanom, family = binomial))
m3 <- with(data, glm(FM ~ Gen_Location+Host+PrecipJlanom, family = binomial))
m4 <- with(data, glm(FM ~ Gen_Location+Host+PrecipAanom, family = binomial))
m5 <- with(data, glm(FM ~ Gen_Location+Host+PrecipJnanom+Year_Collected, family = binomial))
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
AIC(m5)


library(multcomp)
summary(glht(m2, mcp(Gen_Location="Tukey")))
Anova(m2,3)

#to build a general additive model
library(mgcv)
gam1<-gam(FM~s(Day)+Gen_Location+Year_Collected+PrecipJnanom, family = binomial, data=data)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)

gam2<-gam(FM~s(Day)+Gen_Location+Year_Collected+Year_Collected, family = binomial, data=data)
summary(gam2)
plot(gam2)

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

fig<-ggplot(data,aes(x = PrecipJnanom, y = FM))+
  ylab("Proporation female")+
  xlab("May Precip Anomaly")+
  theme_bw()

fig+stat_summary(fun.data = "mean_cl_boot")

#summarizing data
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(plyr)
library(car)
lm1<-glm.nb(Count~Year,dataL)
summary(dataN)
ddply(dataL,"Loc",summarise,Mean=mean(Count),MIN=min(Count), MAX=max(Count),N=length(Count))

confint

#Graphing seasonal variables

m <- glm.nb(Count ~ PrecipAanom, data = dataN)
m
Anova(m)
summary(m)

library(ggplot2)
p<-ggplot(dataA,aes(x = PrecipManom, y = Count))+
  geom_point(aes(),size= 3)+
  ylab("Count")+
  xlab("May Precip Anomaly")+
  scale_y_log10()+
  theme_bw()+
  #geom_abline(intercept=-1.170965,slope=0.040394)+
  geom_abline(intercept=log(1.470965),slope=0.03411)
p

q<-ggplot(dataN,aes(x = PrecipAanom, y = Count))+
  geom_point(aes(color = Host, shape=Loc),size= 3)+
  ylab("Count")+
  xlab("August Precip Anomaly")+
  scale_y_log10()+
  theme_bw()+
  geom_abline(intercept=log(0.6039),
    slope=(0.06734)
  )
q

q<-ggplot(dataN,aes(x = PrecipAanom, y = Count))+
  geom_point(aes(),size= 3)+
  ylab("Count")+
  xlab("August Precip Anomaly")+
  scale_y_log10()+
  theme_bw()+
  geom_abline(intercept=log(0.6039),
              slope=(0.06734)
  )
q

m5 <- with(data, glm(FM ~ Gen_Location+Host+PrecipJnanom, family = binomial))
summary(m5)
m5
fig<-ggplot(data,aes(x = PrecipJnanom, y = FM))+
  ylab("Proporation female")+
  xlab("June Precip Anomaly")+
  theme_bw()+
  geom_abline(intercept=0.3, slope=((0.008336)))
  
#library(Hmisc)
fig+stat_summary(fun.data = "mean_cl_boot")

summary(data)
