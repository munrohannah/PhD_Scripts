Borrelia16sqPCR

setwd("D:/Rfiles/TickPrev") #how to set my working directory
data<-read.csv("Borrelia16sqPCR.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles



#graphing time
library(ggplot2)

p<-ggplot(data,aes(X16s,qpcr))+
  geom_point()+
  theme_minimal()+
  ylim(0,40)
p

p+ylim(0,1000000000000)+theme_minimal()

#correlation
with(data,cor(X16s,qpcr,use="complete.obs"))
with(data,cor.test(X16s,qpcr,use="complete.obs"))

a<-lm(X16s~qpcr,data=data)
a
p+geom_smooth(method = "lm", se = FALSE)
