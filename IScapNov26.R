require(SimRAD)
require(seqinr)

setwd("D:/Rfiles/RsimRAD")


aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)



width(aagen)



ratio <- 1101176065/sum(width(aagen))

ratio





# RAD EcoRI

res1 <- c()

for(i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)
  
  aagen.dig <- insilico.digest(aagen, cut_site_5prime1="G", cut_site_3prime1="AATTC", verbose=FALSE)
  
  ratio <- 1101176065/sum(width(aagen))
  
  cat(ratio, "\n")
  
  res1 <- c(res1, round(length(aagen.dig)*ratio)*2)
  
}

mean(res1)

sd(res1)


## GBS PstI & MspI and ddRAD with size selected 210-260.



#Restriction Enzyme 1

#PstI

cut_site_5prime1 <- "CTGCA"

cut_site_3prime1 <- "G"

#Restriction Enzyme 2

#MspI

cut_site_5prime2 <- "C"

cut_site_3prime2 <- "CGG"



res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 1101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)

