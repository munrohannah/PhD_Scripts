rm(list=ls())

require(SimRAD)
require(seqinr)

setwd("D:/Rfiles/RsimRAD")

rfsq <- ref.DNAseq("D:/Rfiles/RsimRAD/Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fasta", subselect.contigs = TRUE, prop.contigs = 0.10)

# length of the reference sequence:
width(rfsq)
# ratio for the cross-multiplication of the number of fragments and loci at the genomes scale:
genome.size <- 2100000000 # genome size: 2.1Gb
ratio <- genome.size/width(rfsq)
ratio

# computing GC content:
require(seqinr)
GC(s2c(rfsq))

#Restriction Enzyme 1
#PstI
cs_5p1 <- "CTGCA"
cs_3p1 <- "G"

#EcoRI
cs_5p1 <- "G"
cs_3p1 <- "AATTC"


#Restriction Enzyme 2
#MspI #
cs_5p2 <- "C"
cs_3p2 <- "CGG"


simseq.dig <- insilico.digest(rfsq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
simseq.sel <- adapt.select(simseq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)

# wide size selection (200-270):
wid.simseq <- size.select(simseq.sel, min.size = 250, max.size = 300, graph=TRUE, verbose=TRUE)


#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseq)*ratio)
wid.size

# narrow size selection (210-260):
nar.simseq <- size.select(simseq.sel, min.size = 210, max.size = 260, graph=TRUE, verbose=TRUE)

#correcting for the fact that we only looked at a subset of the genome
nar.size <- round(length(nar.simseq)*ratio)

#the resulting fragment characteristics can be further examined:
boxplot(list(width(simseq.sel), width(wid.simseq), width(nar.simseq)), names=c("All fragments",
                                                                               "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")

