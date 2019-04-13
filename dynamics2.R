setwd("D:/Rfiles/TickPrev") #how to set my working directory
Ldata<-read.csv("Larva.csv",header=TRUE) #how to import a file, csv file, and has a header
head(Adata) # gives the first few lines of data
str(Adata) #check the data type and varible
summary(Ldata) #gives means etc

Ndata<-read.csv("Nymphs.csv",header=TRUE)
Adata<-read.csv("Adults.csv",header=TRUE)

#setting up each factor so not treated as a continous variable
Ndata$Year<-factor(Ndata$Year)

library(MASS)
library(car)

m1 <- glm.nb(Count ~ Gen_Location + Year +Host, data = Adata)
m2 <- glm.nb(Count ~ Gen_Location, data = Adata)
m3 <- glm.nb(Count ~ Year, data = Adata)
m4 <- glm.nb(Count ~ Host+Year, data = Ldata)
summary(m1)
Anova(m4,type="III")
anova(m4)

#to get post hoc comparisons
library(multcomp)
summary(glht(m4, mcp(Year="Tukey")))

#to check to see if nb has to be applied or if we can go poisson
var(Adata$Count)/mean(Adata$Count)
hist(Adata$Count, breaks = 50) #histogram

#plot data
plot(Ndata$Day,log10(Ndata$Count))

#to build a general additive model
library(mgcv)
gam1<-gam(Count~s(Day)+Year+Host, family = negbin(theta= c(1,10)), data=Ndata)
summary(gam1)
plot(gam1)
plot(gam1,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(gam1,pages=1,all.terms=TRUE)

gam2<-gam(Count~s(Day), family = negbin(theta= c(1,10)), data=Ldata)
plot(gam2)

#in order to graph the predicted values need to pull out of gam
MD1<-data.frame(Day = seq(from=min(Ldata$Day), to = max(Ldata$Day)))
P1<-(predict(gam2,newdata=MD1,se = TRUE))
SE.UP1<-exp(P1$fit+2*P1$se.fit)
SE.DOWN1<-exp(P1$fit-2*P1$se.fit)
Fit1<-exp(P1$fit)
summary(Fit1)

#this was to predict values on a multi variate data set
library(reshape)
df1 <- data.frame(Host=c("ATPU","COMU","RAZO"))
df2 <- data.frame(Day= min(Ndata$Day):max(Ndata$Day))
df3 <- data.frame(Year= c("2011","2012","2013","2014","2015"))
MD2<-expand.grid.df(df1, df2, df3)

P2<-(predict(gam1,newdata=MD2,se = TRUE, type = "terms"))
P2<-(predict(gam1,newdata=MD2,se = TRUE))
P3<-data.frame(P2)
P4 <- P3[ which(P3$se.fit.Year>1.4 &
                  P3$se.fit.Host >0.7), ]
P4$Day<-(df2$Day)
P4$Fit<-exp(P4$fit.s.Day.)
P4$se.up<-exp(P4$fit.s.Day.+2*P4$se.fit.s.Day.)
P4$se.down<-exp(P4$fit.s.Day.-2*P4$se.fit.s.Day.)

plot(P4$Day,P4$Fit)
lines(P4$Day,P4$se.down)

plot(P3$se.fit.Year)
plot(P3$se.fit.Host)


MD2$SE.UP2<-exp(P2$fit+2*P2$se.fit)
MD2$SE.DOWN2<-exp(P2$fit-2*P2$se.fit)
MD2$Fit2<-exp(P2$fit)

plot(MD2$Day,log10(MD2$Fit2))

library(plyr)
sumdata2 <- ddply(MD2, c("Day") , summarise, #just need to change the factor
                 Fitmean    = mean(Fit2),
                 SEupmean = mean(SE.UP2),
                 SEdownmean = mean(SE.DOWN2))

plot(sumdata2$Fitmean)

sumdata2$Fitmean<-(10^(sumdata2$Fitmean))
sumdata2$SEupmean<-10^(sumdata2$SEupmean)
sumdata2$SEdownmean<-10^(sumdata2$SEdownmean)


plot(sumdata2$Day,sumdata2$Fitmean)
lines(sumdata2$Day,sumdata2$SEupmean)
lines(sumdata2$Day,sumdata2$SEdownmean)

plot(sumdata$Day,sumdata$Fitmean)

#trying something else for graphing gam1
library(visreg)
plotdata <- visreg(gam1, type = "contrast", plot = T)
plotdata2<-data.frame(plotdata)
plot(gam1)


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
  geom_line(aes(x, lower2), linetype = "dashed", P5)+
  geom_line(aes(x, upper2),linetype = "dashed", P5) +
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10(breaks = c(0.01,1,100))+
  scale_x_continuous(breaks = c(151,181,212))+
  theme_bw()

ggplot(MD2,aes(x = Day, y = Fit2))+
  geom_point(aes(color = Host, shape=Year),size= 3)+
  scale_y_log10()+
  ylab("Count")+
  xlab("Julian Day")+
  scale_y_log10()+
  theme_classic()

MDCOMU13<-subset(MD2,Host=="COMU" & Year=="2013")
plot(MDCOMU13$Day,MDCOMU13$Fit2)


library(plyr)
sumdata <- ddply(MD2, c("Day") , summarise, #just need to change the factor
                 Fitmean    = mean(Lmean),
                 SEupmean = mean(SE.UP2),
                 SEdownmean = mean(SE.DOWN2)
)
sumdata

#ggplot
library(ggplot2)
p<-ggplot(Ldata, aes(Day, Count)) + 
  geom_point(aes(color = Host,shape = Year) )+
  theme_classic()
p

p+ geom_line(aes(Day, Fit2), MDCOMU13)+ 
  geom_line(aes(Day, SE.UP2), linetype = "dashed", MDCOMU13)+
  geom_line(aes(Day, SE.DOWN2),linetype = "dashed", MDCOMU13) +
  scale_y_log10()


plot(MD1$Day,(Fit1))
lines(MD2$Day,log10(SE.UP),lty=2)
lines(MD2$Day,log10(SE.DOWN),lty=2)


#now to plot it
plot(Ldata$Day,Ldata$Count)
lines(sumdata$Day,10^(sumdata$Fitmean))
lines(P2[, Day] + attr(P2, "constant"))

lines(MD1$Day,(Fit1))
lines(MD1$Day,log10(SE.UP1),lty=2)
lines(MD1$Day,log10(SE.DOWN1),lty=2)

lines(sumdata$Day,(sumdata$Fitmean))
lines(MD2$Day,log10(SE.UP2),lty=2)
lines(MD2$Day,log10(SE.DOWN2),lty=2)

#now to plot with symbols and colours
plot(Ldata$Day,log10(Ldata$Count),pch=as.integer(Ldata$Year),col=as.integer(Ldata$Host))
lines(MD1$Day,log10(Fit1))
lines(MD1$Day,log10(SE.UP),lty=2)
lines(MD1$Day,log10(SE.DOWN),lty=2)



#box plots for the catagorical
boxplot(log10(Count)~Year,data=Ldata)

#plots of 95CI with box plots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(ggplot2)
plot1 <- ggplot(aes(y = (Count), x = factor(Year)), data = Ldata)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ylab("Count")+
  xlab("Year")+
  scale_y_log10()+
  theme_bw()
plot1

sum(Ndata$Count)
