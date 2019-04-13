rm(list=ls())

require(SimRAD)
require(seqinr)

setwd("D:/Rfiles/RsimRAD")


aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)



width(aagen)



ratio <- 2101176065/sum(width(aagen))

ratio

#Restriction Enzyme 1

#PstI

cut_site_5prime1 <- "CTGCA"

cut_site_3prime1 <- "G"

#Restriction Enzyme 2

#MspI

cut_site_5prime2 <- "C"

cut_site_3prime2 <- "CGG"

#Restriction Enzyme 3

#EcoRI

cut_site_5prime3 <- "G"

cut_site_3prime3 <- "AATTC"

#Restriction Enzyme 4

#MseI

cut_site_5prime4 <- "T"

cut_site_3prime4 <- "TAA"


## GBS PstI & MspI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 2101176065/sum(width(aagen))
  
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


## GBS PstI & MspI and ddRAD with size selected 110-160.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 2101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 110, max.size = 160, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)


## GBS PstI & MspI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 2101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime3, cut_site_3prime3, cut_site_5prime2, cut_site_3prime2, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)

## GBS MspI & EcoRI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 2101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime2, cut_site_3prime2, cut_site_5prime3, cut_site_3prime3, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)

## GBS PstI & EcoRI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 2101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime3, cut_site_3prime3, cut_site_5prime1, cut_site_3prime1, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)


## GBS MseI & MspI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 1101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime4, cut_site_3prime4, cut_site_5prime2, cut_site_3prime2, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)

## GBS MseI & PstI and ddRAD with size selected 210-260.

res4 <- c()

res5 <-c()

for (i in 1:30)
  
{
  
  cat(i,"; ")
  
  aagen <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)    
  
  ratio <- 1101176065/sum(width(aagen))
  
  aagen.select <- adapt.select(insilico.digest(aagen, cut_site_5prime4, cut_site_3prime4, cut_site_5prime1, cut_site_3prime1, verbose=FALSE), type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2)
  
  aagen.sized <- size.select(aagen.select, min.size = 210, max.size = 260, graph=FALSE, verbose=FALSE)
  
  cat(ratio, "\n")
  
  res4 <- c(res4, round(length(aagen.select)*ratio))
  
  res5 <- c(res5, round(length(aagen.sized)*ratio))
  
}

mean(res4)

sd(res4)

mean(res5)

sd(res5)


boxplot(list(width(simseq.sel), width(wid.simseq), width(nar.simseq)), names=c("All fragments",
                                                                               "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")
