rm(list=ls())

require(SimRAD)
require(seqinr)

setwd("D:/Rfiles/RsimRAD")

rfsq <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_CONTIGS_IscaW1.fa", subselect.contigs = TRUE, prop.contigs = 0.10)

# length of the reference sequence:
width(rfsq)
# ratio for the cross-multiplication of the number of fragments and loci at the genomes scale:
genome.size <- 2100000000 # genome size: 2.1Gb
ratio <- genome.size/width(rfsq)
ratio

#Restriction Enzyme
#PstI
cs_5pPstI <- "CTGCA"
cs_3pPstI <- "G"

#EcoRI
cs_5pEcoRI <- "G"
cs_3pEcoRI <- "AATTC"


#MspI
cs_5pMspI <- "C"
cs_3pMspI <- "CGG"

#MluCl
cs_5pMluCl <- "AATT"
cs_3pMluCl <- ""

#SphI
cs_5pSphI <- "GCATG"
cs_3pSphI <- "C"

#NlaIII
cs_5pNlaIII <- "CATG"
cs_3pNlaIII<- ""

#SphI and EcoRI digest
simseqSphI.EcoRI.dig <- insilico.digest(rfsq, cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI, verbose=TRUE)
simseqSphI.EcoRI.sel <- adapt.select(simseqSphI.EcoRI.dig, type="AB+BA", cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI)
# wide size selection (200-270):
wid.simseqSphI.EcoRI <- size.select(simseqSphI.EcoRI.sel, min.size = 300, max.size = 350, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqSphI.EcoRI)*ratio)
wid.size
