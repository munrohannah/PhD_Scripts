#rarifaction curve for MLST ST and alleles

library(vegan)
library(tidyverse)
library(dplyr)
library(tidyr)

setwd("D:/Rfiles/MLST") 
data<-read.csv("ST.csv",header=TRUE) #this is a file with just ST
data<-select(data,ST244:ST694) #removing the rows with id
names(data)

##turns out what i want is a species accumlation curve
##code from https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/specaccum

sp1 <- specaccum(data)
sp1
sp21 <- specaccum(data, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")
plot(sp2)

#loading the allele data
dataA<-read.csv("mlstALL3.csv",header=TRUE) #this is a file allele numbers by id
#cleaning data for use and some wrangling
data2<-gather(dataA,"individual",2:9)
data2$value<-1
summary (data2)
names(data2)[1]<-"ID"
names(data2)[3]<-"alleleNO"
library(reshape2)
data3<-dcast(data2,ID~alleleNO,mean)
summary(data3)
head(data3)
data3[data3 == "NaN"] <- 0
data3<-select(data3,clpa116:uvra39)

sp1 <- specaccum(data3)
sp1
sp22 <- specaccum(data3, "random")
sp22
summary(sp22)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp22, col="yellow", add=TRUE, pch="+")

##Trying something new to graphs this is what I used
# Combine the specaccum objects into a list 
l <- list(sp21, sp22) 

# Calculate required y-axis limits
ylm <- range(sapply(l, '[[', 'richness') + 
               sapply(l, '[[', 'sd') * c(-2, 2))

# Apply a plotting function over the indices of the list
sapply(seq_along(l), function(i) {
  if (i==1) { # If it's the first list element, use plot()
    with(l[[i]], {
      plot(sites, richness, type='l', ylim=ylm, 
           xlab='Samples', ylab='Richness', las=1)
      segments(seq_len(max(sites)), y0=richness - 1*sd, 
               y1=richness + 1*sd)
    })    
  } else {
    with(l[[i]], { # for subsequent elements, use lines()
      lines(sites, richness, col=i)
      segments(seq_len(max(sites)), y0=richness - 1*sd, 
               y1=richness + 1*sd, col=i)
    })     
  }
})

legend('topleft', c('STs', 'Alleles'), col=1:2, lty=1, 
       bty='n'
       , inset=0.005)
