#entering female data with averages in missing
setwd("D:/Rfiles/ticks") #how to set my working directory
##measure_2017_com is female ##male is male data
data<-read.csv("measure_2017_com.csv",header=TRUE)#how to import a file, csv file, and has a header
summary(data)

#female is 11:34 #male is 10:36
v<-data[,11:34]
library("Hmisc")
res2 <- rcorr(as.matrix(log(v)))
res2
res <- cor(v)
library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#log transform the data
#female is 11:34 #male is 10:36
datalog<-log(data[,11:34])

#get the shapevalues
#female is 11:34 #male is 10:36
vars <- names(data)[11:34]
fits <- lapply(vars, function(x) {glm(substitute(i ~ Host, list(i = as.name(x))), data = data)})
lapply(fits,summary) # this works
lapply(fits, coefficients) # this works
lapply(fits, residuals)

datalogS<-data.frame(lapply(fits, residuals))

#PCA ##change the data file name ###datalog or datalogS
res.pca <- prcomp(datalog, scale = TRUE)
library(factoextra)
var <- get_pca_var(res.pca)
var
#contribution of each dimention
k<-data.frame(var$contrib[, 1:4])
k

#adding the PC to file
scores<-data.frame(res.pca$x)
data$PC1log<-scores$PC1
data$PC2log<-scores$PC2
data$PC3logS<-scores$PC3

#contributions #with lables removed
fviz_pca_biplot(res.pca, col.var="contrib",
                label= "none"
                )+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=3) + theme_minimal()


#AIC
m1<-glm(PC1log~Gen_Location+Year+Host, data= data)
m2<-glm(PC1log~Year+Host, data= data)
m3<-glm(PC1log~Host+Gen_Location, data= data)
m4<-glm(PC1log~Gen_Location+Year, data= data)
m5<-glm(PC1log~Year, data= data)
m6<-glm(PC1log~Host, data= data)
m7<-glm(PC1log~Gen_Location, data= data)
AIC(m1,m2,m3,m4,m5,m6,m7)

m10<-glm(PC2log~Gen_Location+Host+Pos_Neg, data=data)
AIC(m10)
logLik(m10)

##GAM time
library(mgcv)
gam1<-gam(PC2logS~s(Jday)+Host+Gen_Location, data=data)
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
plot1 <- ggplot(aes(y = (PC2log), x = factor(Gen_Location)), data = data)
plot1 <- plot1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3,aes(color = factor(Host))) + 
  ylab("PC1")+
  xlab("Year")+
  theme_bw()
plot1

