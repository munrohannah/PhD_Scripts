### with multple variables
setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("Measure2_rAnaysisSept14.csv",header=TRUE)#how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles
library(reshape2)

cdata<-dcast(data,Tick_Vial+ Age+Engorged+Gen_Location+Date_Collected+Host+Pos_Neg+Measured_by~Measurement,mean)
head(cdata)

flat<-cdata[cdata$Engorged=="N",]
head(flat)

female<-flat[flat$Age=="AF",]
male<-flat[flat$Age=="AM",]
summary(male)

write.csv(female, file = "femaleALL.csv")
write.csv(male, file = "maleALL.csv")

#####change depending 
setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("maleALL.csv",header=TRUE)#how to import a file, csv file, and has a header
summary(data)
data<-na.exclude(data)



#female is 11:34 #male is 10:36
v<-data[,11:14]
library("Hmisc")
res2 <- rcorr(as.matrix(log(v)))
res2
res <- cor(v)
library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(res)

#log transform the data
datalog<-log(data[,11:14])

#get the shapevalues
vars <- names(data)[11:14]
fits <- lapply(vars, function(x) {glm(substitute(i ~ Host, list(i = as.name(x))), data = data)})
lapply(fits,summary) # this works
lapply(fits, coefficients) # this works
lapply(fits, residuals)

datalogS<-data.frame(lapply(fits, residuals))

#PCA ##change the data file name ###datalog or datalogS
res.pca <- prcomp(datalogS, scale = TRUE)
library(factoextra)
var <- get_pca_var(res.pca)
var
#contribution of each dimention
k<-data.frame(var$contrib[, 1:4])
k

#adding the PC to file
scores<-data.frame(res.pca$x)
data$PC1logS<-scores$PC1
data$PC2logS<-scores$PC2
data$PC3log<-scores$PC3
data$PC4log<-scores$PC4
#contributions #with lables removed
fviz_pca_biplot(res.pca, col.var="contrib",
                label= "var"
)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=24) + theme_minimal()
data$Year<-as.factor(data$Year)

#AIC
m1<-glm(PC1logS~Gen_Location+Year+Host+Measured_by, data= data)
m2<-glm(PC1logS~Year+Host+Measured_by, data= data)
m3<-glm(PC1logS~Host+Gen_Location+Measured_by, data= data)
m4<-glm(PC1logS~Gen_Location+Year+Measured_by, data= data)
m5<-glm(PC1logS~Year+Measured_by, data= data)
m6<-glm(PC1logS~Host+Measured_by, data= data)
m7<-glm(PC1logS~Gen_Location+Measured_by, data= data)
m8<-glm(PC1logS~Measured_by, data= data)
AIC(m1,m2,m3,m4,m5,m6,m7,m8)
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)
logLik(m6)
logLik(m7)
logLik(m8)

m10<-glm(PC1logS~Year+Measured_by+Pos_Neg, data=data)
AIC(m10)
logLik(m10)

m11<-glm(PC1logS~Year, data=data)
AIC(m11)
logLik(m11)

#if i remove measured-by
m1<-glm(PC1log~Year, data= data)
m2<-glm(PC2log~Gen_Location, data= data)
m3<-glm(PC1logS~Year, data= data)
m4<-glm(PC2logS~Year, data= data)
AIC(m1,m2,m3,m4)

library(car)
Anova(m4)

##GAM time
library(mgcv)
gam1<-gam(PC2logS~s(Jday)+Year, data=data)
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

P5 <- smooths[ which(smooths$Variable== "Jday"), ]

plot(P5$x,P5$smooth)
lines(P5$x,P5$lower)
lines(P5$x,P5$upper)

#graphing with lines
library(ggplot2)
p<-ggplot(data,aes(x = Jday, y = PC2logS))+
  geom_point(aes(color = Host, shape= Gen_Location) ,size= 3)+
  geom_line(aes(x, smooth), P5)+
  ylab("PC1")+
  xlab("Julian Day")+
  #scale_y_log10(limit=c(0.15,120))+
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
plot1 <- ggplot(aes(y = (PC2logS), x = factor(Year)), data = data)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3,aes(color = factor(Host))) + 
  ylab("PC1")+
  xlab("Year")+
  theme_bw()
plot1
