#borrelia prev in tick

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("Nov3data.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
data$Year<-factor(data$Year)
library(MASS)
library(car)
m1 <- with(data, glm(Pos_Neg ~ Host, family = binomial))
summary(m1)
Anova(m1,3)


#serology

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("ser.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

#setting up each factor so not treated as a continous variable
data$pos<-factor(data$pos)
data$Year<-factor(data$Year)

#GLM time
library(MASS)
library(car)
m1 <- with(data, glm(pos ~ Species, family = binomial))
m2 <- with(data, glm(pos ~ Loc, family = binomial))
m3 <- with(data, glm(pos ~ Year, family = binomial))
m4 <- with(data, glm(pos ~ Species+Loc, family = binomial))
m5 <- with(data, glm(pos ~ Species+Year, family = binomial))
m6 <- with(data, glm(pos ~ Loc+Year, family = binomial))
m7 <- with(data, glm(pos ~ Loc+Year+Species, family = binomial))
AIC(m1,m2,m3,m4,m5,m6,m7)
summary(m1)
Anova(m2,3)


library(multcomp)
summary(glht(m1, mcp(Species="Tukey")))

library(ggplot2)
fig2<-ggplot(data,aes(x=Species, y=pos))+
  stat_summary(fun.data = 'mean_cl_boot')+
  ylab("seroprevalence")+
  xlab("species")+
  scale_y_continuous(limit=c(0,1),breaks = c(0,0.5,1))+
  theme_bw()
fig2

#gam to see the influence of date
library(mgcv)
gam1<-gam(pos~s(Jdate)+Species, family = binomial, data=data)
summary(gam1)
plot(gam1)

fig2+stat_summary(fun.data = "mean_cl_boot",geom = "boxplot")
fig+stat_smooth()