library(SNPassoc)

setwd("D:/Rfiles/Ticks") #how to set my working directory
data<-read.csv("batch10c.csv",header=TRUE) #how to import a file, csv file, and has a header
x<-summary(data) #gives means etc
x

#put all the SNPs togeather

myData<-setupSNP(data=data,colSNPs=5:222,sep="/")
summary(myData)
snp<-data[,5:222]

plotMissing(myData)
res<-tableHWE(myData)
res<-data.frame(res)
res<- tableHWE(myData,strata=myData$Island)
res

association(run~X266,data=myData)

ans<-WGassociation(Island~1+run,data=myData,model="all")
summary(ans)
ans
ans$codominant
res$codominate<-ans$codominant

write.csv(res, file = "hwe_etc.csv", row.names=T)
write.csv(res, file = "hwe.csv")

plot(ans,cex=0.8)
