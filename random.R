emptyvar<-numeric(1000)                     #build an empty variable for the prevalence to be placed in
for (i in 1:1000) {                         #starting a loop
  emptyvar[i]<-mean(sample(AE2,replace=T),) #place mean of randomly sampled data, must replace vector
}                                           #close the loop
AE2CI<-quantile(emptyvar, c(0.025, 0.975))   #give the 95%CI

allele<-data.frame(c("a","b","a","b","c","a","c","b","a","c","c","a","b","d"))
rare<-numeric(100)
for (i in 1:100) {
  rare[i]<-dim(unique(sample(allele,5,replace=T),))
}
