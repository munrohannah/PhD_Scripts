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

#MluCl and EcoRI digest
simseqMluCl.EcoRI.dig <- insilico.digest(rfsq,  cs_5pEcoRI, cs_3pEcoRI,cs_5pMluCl, cs_3pMluCl, verbose=TRUE)
simseqMluCl.EcoRI.sel <- adapt.select(simseqMluCl.EcoRI.dig, type="AB+BA",  cs_5pEcoRI, cs_3pEcoRI,cs_5pMluCl, cs_3pMluCl)
# wide size selection (200-270):
wid.simseq <- size.select(simseqMluCl.EcoRI.sel, min.size = 100, max.size = 150, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseq)*ratio)
wid.size

#NlaIII and EcoRI digest
simseqNlaIII.EcoRI.dig <- insilico.digest(rfsq, cs_5pEcoRI, cs_3pEcoRI,cs_5pNlaIII, cs_3pNlaIII, verbose=TRUE)
simseqNlaIII.EcoRI.sel <- adapt.select(simseqNlaIII.EcoRI.dig, type="AB+BA",  cs_5pEcoRI, cs_3pEcoRI,cs_5pNlaIII, cs_3pNlaIII)
# wide size selection (200-270):
wid.simseq <- size.select(simseqNlaIII.EcoRI.sel, min.size = 300, max.size = 350, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseq)*ratio)
wid.size

#MluCl and NlaIII digest
simseqMluCl.NlaIII.dig <- insilico.digest(rfsq, cs_5pMluCl, cs_3pMluCl, cs_5pNlaIII, cs_3pNlaIII, verbose=TRUE)
simseqMluCl.NlaIII.sel <- adapt.select(simseqMluCl.NlaIII.dig, type="AB+BA", cs_5pMluCl, cs_3pMluCl, cs_5pNlaIII, cs_3pNlaIII)
# size selection (200-270):
wid.simseqMluCl.NlaIII <- size.select(simseqMluCl.NlaIII.sel, min.size = 300, max.size = 350, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqMluCl.NlaIII)*ratio)
wid.size

#MluCl and SphI digest
simseqMluCl.SphI.dig <- insilico.digest(rfsq, cs_5pMluCl, cs_3pMluCl, cs_5pSphI, cs_3pSphI, verbose=TRUE)
simseqMluCl.SphI.sel <- adapt.select(simseqMluCl.SphI.dig, type="AB+BA", cs_5pMluCl, cs_3pMluCl, cs_5pSphI, cs_3pSphI)
# size selection (200-270):
wid.simseqMluCl.SphI <- size.select(simseqMluCl.SphI.sel, min.size = 100, max.size = 150, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqMluCl.SphI)*ratio)
wid.size

#SphI and NlaIII digest
simseqSphI.NlaIII.dig <- insilico.digest(rfsq, cs_5pSphI, cs_3pSphI, cs_5pNlaIII, cs_3pNlaIII, verbose=TRUE)
simseqSphI.NlaIII.sel <- adapt.select(simseqSphI.NlaIII.dig, type="AB+BA", cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI)
# wide size selection (200-270):
wid.simseqSphI.NlaIII <- size.select(simseqSphI.NlaIII.sel, min.size = 200, max.size = 250, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqSphI.NlaIII)*ratio)
wid.size

#SphI and EcoRI digest
simseqSphI.EcoRI.dig <- insilico.digest(rfsq, cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI, verbose=TRUE)
simseqSphI.EcoRI.sel <- adapt.select(simseqSphI.EcoRI.dig, type="AB+BA", cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI)
# wide size selection (200-270):
wid.simseqSphI.EcoRI <- size.select(simseqSphI.EcoRI.sel, min.size = 300, max.size = 350, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqSphI.EcoRI)*ratio)
wid.size

#PstI and MspI digest
simseqPstI.MspI.dig <- insilico.digest(rfsq, cs_5pPstI, cs_3pPstI, cs_5pMspI, cs_3pMspI, verbose=TRUE)
simseqPstI.MspI.sel <- adapt.select(simseqPstI.MspI.dig, type="AB+BA", cs_5pPstI, cs_3pPstI, cs_5pEcoRI, cs_3pEcoRI)
# wide size selection (200-270):
wid.simseqPstI.MspI <- size.select(simseqPstI.MspI.sel, min.size = 250, max.size = 300, graph=TRUE, verbose=TRUE)
#correcting for the fact that we only looked at a subset of the genome
wid.size <- round(length(wid.simseqPstI.MspI)*ratio)
wid.size


#
simseqPstI.MspI.dig <- insilico.digest(rfsq, cs_5pPstI, cs_3pPstI, cs_5pMspI, cs_3pMspI, verbose=TRUE)
simseqSphI.EcoRI.dig <- insilico.digest(rfsq, cs_5pSphI, cs_3pSphI, cs_5pEcoRI, cs_3pEcoRI, verbose=TRUE)
simseqSphI.NlaIII.dig <- insilico.digest(rfsq, cs_5pSphI, cs_3pSphI, cs_5pNlaIII, cs_3pNlaIII, verbose=TRUE)
simseqMluCl.NlaIII.dig <- insilico.digest(rfsq, cs_5pMluCl, cs_3pMluCl, cs_5pNlaIII, cs_3pNlaIII, verbose=TRUE)
simseqMluCl.EcoRI.dig <- insilico.digest(rfsq,  cs_5pEcoRI, cs_3pEcoRI,cs_5pMluCl, cs_3pMluCl, verbose=TRUE)
