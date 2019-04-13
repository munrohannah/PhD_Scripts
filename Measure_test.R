measure_test

setwd("D:/Rfiles/ticks") #how to set my working directory
data<-read.csv("measure_2017.csv",header=TRUE)#how to import a file, csv file, and has a header
com<-read.csv("measure_2017_com.csv",header=TRUE)
comM<-read.csv("measure_2017_comM.csv",header=TRUE)
com$BCL<-as.numeric(com$BCL)
com$year<-as.factor(com$year)
hist(comM$C1W)

scu<-read.csv("scutal.csv",header=TRUE)
scuM <- scu[ which(scu$Age=='AM'), ]
scuF <- scu[ which(scu$Age=='AF'), ]

ot <- scuM[ which(scuM$Measured.by=='ot'), ]
scutal<-with(ot,cbind(SCL,SCW))


library(MAINT.Data)
m1<-manova(scutal~Gen_Location+Host, data= ot)

m2<-MANOVA(scutal~Host, data= ot)
m3<-MANOVA(scutal~Gen_Location, data= ot)
m4<-MANOVA(scutal~qPCR_Result, data= ot)
m5<-MANOVA(scutal~Measured.by, data= ot)
m6<-MANOVA(scutal~Measured.by+qPCR_Result, data= ot)
m7<-MANOVA(scutal~Engorged, data= ot)
m8<-MANOVA(scutal~Host+Engorged, data= ot)
m9<-MANOVA(scutal~Host+Engorged+Measured.by, data= ot)

AIC(m1)
library(car)
Anova(m4)

comU <- com[ which(com$Host!='Unknown'), ]

summary(com)
comVAR<-comU[,11:35]

comM<-comM[1:20,1:35]
comVAR<-comM[,9:35]

comVARlog<-log(comVAR[,])


means<-aggregate(com[, 11:35], list(com$Host), mean)
means

a<-glm(BCL~Host,data=com)
var.names <- c(com[,11:35])

vars <- names(com)[11:35]
fits <- lapply(vars, function(x) {glm(substitute(i ~ Host, list(i = as.name(x))), data = com)})
lapply(fits,summary) # this works
lapply(fits, coefficients) # this works
a<-data.frame(lapply(fits, residuals))
#lapply(fits, summary(fits)$coefficients[,4])# this for example does not work


res.pca <- prcomp(comVAR, scale = TRUE)
names(res.pca)
b<-data.frame(res.pca$rotation)
res.pca$scale
res.pca$center
head(res.pca$sdev)
head(unclass(res.pca$rotation)[, 1:4])
head(unclass(res.pca$rotation))
# Eigenvalues
eig <- (res.pca$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.com <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.com)

barplot(eig.com[, 2], names.arg=1:nrow(eig.com), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.com), 
      eig.com[, 2], 
      type="b", pch=19, col = "red")
library(factoextra)
var <- get_pca_var(res.pca)
var
var$coord[, 1:4]
var$contrib[, 1:4]


scores<-data.frame(res.pca$x)
comM$PC1<-scores$PC1
comM$PC2<-scores$PC2
comM$PC3<-scores$PC3


fviz_pca_var(res.pca, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=3) + theme_minimal()

fviz_pca_var(res.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=2) + theme_minimal()

ind.coord <- res.pca$x
head(ind.coord[, 1:4])
# Compute the square of the distance between an individual and the
# center of gravity
center <- res.pca$center
scale<- res.pca$scale
getdistance <- function(ind_row, center, scale){
  return(sum(((ind_row-center)/scale)^2))
}
d2 <- apply(comVAR,1,getdistance, center, scale)
# Compute the cos2
cos2 <- function(ind.coord, d2){return(ind.coord^2/d2)}
ind.cos2 <- apply(ind.coord, 2, cos2, d2)
head(ind.cos2[, 1:4])

fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.50) + theme_minimal()
fviz_pca_biplot(res.pca,  geom = "text") +
  theme_minimal()

scores<-data.frame(pc.cr$scores)
com$PC2<-scores$Comp.1
com$PC2<-scores$Comp.2
com$PC3<-scores$Comp.3

pc.cr$loadings
loadings(pc.cr)


m1<-glm(PC1a~Gen_Location+year+Host, data= com)
m2<-glm(PC1a~year+Host, data= com)
m3<-glm(PC1a~Host+Gen_Location, data= com)
m4<-glm(PC1a~Gen_Location+year, data= com)
m5<-glm(PC1a~year, data= com)
m6<-glm(PC1a~Host, data= com)
m7<-glm(PC1a~Gen_Location, data= com)
m8<-glm(PC1a~Pos_Neg, data= com)
m9<-glm(PC1a~year+Host+Pos_Neg, data= com)

m1<-glm(PC2~Gen_Location+Host, data= comM)
m2<-glm(PC2~Host, data= comM)
m3<-glm(PC2~Host+Gen_Location, data= comM)
m4<-glm(PC2~Gen_Location, data= comM)
m5<-glm(PC2~year, data= comM)
m6<-glm(PC2~Host, data= comM)
m7<-glm(PC2~Gen_Location, data= comM)
m8<-glm(PC2~qPCR_Result, data= comM)
m9<-glm(PC2~Host+Gen_Location+qPCR_Result, data= comM)

AIC(m1)
AIC(m2)
AIC(m3)
AIC(m4)
AIC(m5)
AIC(m6)
AIC(m7)
AIC(m8)
AIC(m9)

library(car)
Anova(m9)

library(MASS)
library(car)
library(mgcv)
library(visreg)
library(plyr)
library(ggplot2)
library(boot)

gamT<-gam(PC2log~s(Jday)+year+Host+Pos_Neg, data=comU)


summary(gamT)
plot(gamT)
plotdataT <- visreg(gamT, type = "contrast", plot = F)

smoothsT <- ldply(plotdataT, function(part)   
  data.frame(Variable = part$meta$x, 
             x=part$fit[[part$meta$x]], 
             smooth=part$fit$visregFit, 
             lower=part$fit$visregLwr, 
             upper=part$fit$visregUpr))

P5T <- smoothsT[ which(smoothsL$Variable== "Jday"), ]


P5T$smooth2 <- P5T$smooth


p<-ggplot(comU,aes(x = Jday, y = PC2log))+
  geom_point(aes(color = Host, shape= Gen_Location),size= 3)+
  geom_line(aes(x, smooth2), P5T)+
  ylab("Count")+
  xlab("Julian Day")+
  theme_bw()
p

#graphing with lines
library(ggplot2)

fig2<-ggplot(comU,aes(x=Gen_Location, y=PC2log))+
  stat_summary(fun.data = 'mean_cl_boot')+
  ylab("logPC2")+
  xlab("Host")+
  scale_y_continuous()+
  theme_bw()
fig2

scutplot<-ggplot(scu,aes(x=SCL, y=SCW))+
  geom_point(aes(color = Measured.by, shape= Age),size= 3)+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_bw()
scutplot
