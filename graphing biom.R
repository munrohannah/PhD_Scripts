setwd("D:/Rfiles/16s") #how to set my working directory
dataR<-read.csv("L6_new.csv",header=TRUE) #how to import a file, csv file, and has a header
head(dataR) # gives the first few lines of data
str(dataR) #check the data type and varibles

library(reshape2)
mdR <- melt(dataR, id=c("order","bacteria")) #melt the dataso that each observation is seperate

#write to tab
write.csv(mdR, file = "L6_melt.csv",row.names=TRUE)
mdR<-read.csv("L6_melt.csv",header=TRUE) #how to import a file, csv file, and has a header


library(ggplot2)

  
#need to change order of names so that they have neg then pos  
mdR$variable <- factor(mdR$variable, levels = mdR$variable[order(mdR$orderSam)])
mdR$bacteria <- factor(mdR$bacteria, levels = mdR$bacteria[order(mdR$orderOTU)])
  
p<-ggplot(mdR, aes(x = variable, y = value)) +
  geom_bar(aes(fill = bacteria), stat="identity", position="stack") +
  theme_minimal()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
p+scale_fill_manual(values=colours1)
p

colours1<-c("#800000","#1f78b4","#808000","#33a02c","#ffff00",
            "#e31a1c","#008000","#ff7f00","#00ff00","#800080",
            "#1f78b4","#ff00ff","#33a02c","#000080","#e31a1c",
            "#ff6600","#33a02c","#008080","#6a3d9a","#c83737",
            "#b15928","#aad400","#ff7f00")








##L2

setwd("D:/Rfiles/16s") #how to set my working directory
dataR<-read.csv("L2_new.csv",header=TRUE) #how to import a file, csv file, and has a header
head(dataR) # gives the first few lines of data
str(dataR) #check the data type and varibles

library(reshape2)
mdR <- melt(dataR, id=c("order","Bacteria")) #melt the dataso that each observation is seperate

#write to tab
write.csv(mdR, file = "L2_melt.csv",row.names=TRUE)
mdR<-read.csv("L2_melt.csv",header=TRUE) #how to import a file, csv file, and has a header


library(ggplot2)


#need to change order of names so that they have neg then pos  
mdR$variable <- factor(mdR$variable, levels = mdR$variable[order(mdR$orderSam)])
mdR$bacteria <- factor(mdR$bacteria, levels = mdR$bacteria[order(mdR$orderOTU)])

p<-ggplot(mdR, aes(x = variable, y = value)) +
  geom_bar(aes(fill = Bacteria), stat="identity", position="stack") +
  theme_minimal()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p+scale_fill_manual(values=colours1)
p

colours1<-c("#800000","#1f78b4","#808000","#33a02c","#ffff00",
            "#e31a1c","#008000","#ff7f00","#00ff00","#800080",
            "#1f78b4","#ff00ff","#33a02c","#000080","#e31a1c",
            "#ff6600","#33a02c","#008080","#6a3d9a","#c83737",
            "#b15928","#aad400","#ff7f00")

