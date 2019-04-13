setwd("D:/Rfiles/16s") #how to set my working directory
data2<-read.csv("weig_PCoA_ave_joost.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data2) # gives the first few lines of data
str(data2) #check the data type and varibles

library(reshape2)# this is for reshaping collated data
md <- melt(data2, id=(c("sample", "PC"))) #melt the data so that each observation is seperate
md2 <-dcast(md,sample+variable~PC) #cast so that the data has PC1 and PC2 as seperate rows

#graphing collated data
library(ggplot2)
p <- ggplot(md2, aes(PC1, PC2))
p +  geom_point(aes(colour = sample))

#graphing summarized data
library(ggplot2)
p <- ggplot(data, aes(PC1_a, PC2_a, shape = factor(Loc)))
q<-p +  geom_point(aes(size = Bg))
q + geom_point(aes(colour = factor(Age)))

ggplot(data2, aes(x=PC1_a,y=PC2_a))+
  geom_point(aes(color=Age,size=Bg_N,shape=Loc))+ 
  scale_size(breaks=2, range = c(2, 4))+
  theme_minimal()

