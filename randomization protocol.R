#calaculating 95CI through randomization

setwd("D:/Rfiles/16s") #how to set my working directory
dataR<-read.csv("L2_top1per.csv",header=TRUE) #how to import a file, csv file, and has a header
head(data) # gives the first few lines of data
str(data) #check the data type and varibles

#pulling out the pos_neg data based on differnet factors.
Prot<-dataR[,5] ##must change data name, factor name, and specific factor, 14 refers to pos_neg column

#randomization calcluation of 95CI
emptyvar<-numeric(1000)                     #build an empty variable for the prevalence to be placed in
for (i in 1:1000) {                         #starting a loop
  emptyvar[i]<-mean(sample(Prot,replace=T),) #place mean of randomly sampled data, must replace vector
}                                           #close the loop
ProtCI<-quantile(emptyvar, c(0.025, 0.975))   #give the 95%CI


#appending CI to the end of the cdata dataframe
bactCI<-data.frame(t(data.frame(FirCI,SpiCI,ActiCI,BactCI,ProtCI)))
