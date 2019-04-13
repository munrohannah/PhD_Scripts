####Based on a tutorial prepared from 
1#https://botany.natur.cuni.cz/hodnocenidat/Lesson_05_tutorial.pdf

library("vcfR") 
library("adegenet") 
library("adegraphics") 
library("pegas") 
library("StAMPP") 
library("lattice") 
library("gplots") 
library("ape") 
library("ggmap") 

setwd("D:/Rfiles/ddrad")

meta<-read.csv("META.csv",header=TRUE)

vcf <- read.vcfR("NEW3.recode2.vcf")   #read in all data 
#vcf <- read.vcfR("batch_11_b.vcf")   #read in all data 
head(vcf)               #check the vcf object
vcf@fix[1:5,1:8]       #check 



#chrom <- create.chromR(vcf) 
#plot(chrom) # plot the data 

heatmap.bp(vcf)
summary(a$fix)

#quick check read depth distribution per individual 
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
summary(dp)
head(dp)
pdf("DP_RAD_data.pdf", width = 10, height=3) # boxplot 
par(mar=c(8,4,1,1)) 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", 
        las=2, cex=0.4, cex.axis=0.5)
dev.off()

summary(dp)

pdf("DP_RAD_data_zoom.pdf", width = 10, height=3) # boxplot 
par(mar=c(8,4,1,1)) 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", 
        las=2, cex=0.4, cex.axis=0.5, ylim=c(0,50)) 
abline(h=4, col="red") 
dev.off() 

### convert to genlight  
aa.genlight <- vcfR2genlight(vcf, n.cores=1) 
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names 
strata(aa.genlight)<-meta
setPop(aa.genlight)<-~island

# check the genlight 
aa.genlight                        # check the basic info on the genlight object 
indNames(aa.genlight)              # check individual names 
as.matrix(aa.genlight)[1:3,1:3]  # see tiny bit of the data 
pop(aa.genlight)                   # population assignment


# look at the total data matrix (0,1,2; white = missing data) 
glPlot (aa.genlight)  # takes some time 
glPlot (aa.genlight.z)

# N missing SNPs per sample 
x <- summary(t(as.matrix(aa.genlight))) 
x
write.table(x[7,], file = "missing.persample.txt", sep = "\t")  # NAs, if present, are in seventh row of summary 

aa.genlight.z <- new("genlight", (as.matrix(aa.genlight)) 
                     [,(colSums(as.matrix(aa.genlight)) > 0) &  
                         (colSums(is.na(as.matrix(aa.genlight))) == 
                            0)]) # remove the reference-only positions AND remove columns with NA 
aa.genlight.z

###plot total AFS of the dataset 
mySum <- glSum(aa.genlight, alleleAsUnit = TRUE)  
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  
        main="Distribution of ALT allele counts in total dataset") 

##plot AFS per one pop  
aa.genlight.sep <- seppop(aa.genlight, drop=TRUE)  #separate genlights per population 
aa.genlight.sep$GULL

# after seppop you must remove the nonvariant positions within the population 
n.alleles.k <-colSums(as.matrix(aa.genlight.sep$GULL)) # how many alternative alleles are in each locus? 
summary(as.factor(n.alleles.k))                            # how many particular categories of alternative allele counts are in my pop? 
aa.genlight.k <- new("genlight", (as.matrix(aa.genlight.sep$GULL)) 
                       [,(colSums(as.matrix(aa.genlight.sep$GULL)) > 0) &  
                           (colSums(is.na(as.matrix(aa.genlight.sep$GULL))) == 
                              0)]) # remove the reference-only positions AND remove columns with NA 
aa.genlight.k 
summary(colSums(as.matrix(aa.genlight.k)))  #  check if there are no zeros 
# plot unfolded AFS - for one pop. 
mySum <- glSum(aa.genlight.k, alleleAsUnit = TRUE)                     
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  
        main="Distribution of ALT allele counts in BEL")    # plot the original counts of each category 


#### plot AFS for all pops in a batch 
aa.genlight.sep <- seppop(aa.genlight, drop=TRUE)  # separate genlight per population 
# remove the nonvariant positions AND columns with NA within that pop. 
aa.genlight.sep.2 <- lapply (aa.genlight.sep, function (pop) 
{new("genlight", (as.matrix(pop))[,(colSums(as.matrix(pop)) > 0) 
                                  & (colSums(is.na(as.matrix(pop))) == 0)])})  
##add pop identity to list elements 
listnames<-names(aa.genlight.sep.2) 
for (i in seq(listnames)) {pop(aa.genlight.sep.2[[i]])<-
  strata((aa.genlight.sep.2[[i]]),~island)}  
# loop over each population in a list of populations and draw AFS into one fig 
pdf("AFS_all_barplot.pdf", width=5, height=5) 
par(mfrow=c(2,3),mar=c(2,2,2,0)) 
mySum <- lapply (aa.genlight.sep.2, function (pop) { 
  barplot(table(glSum(pop, alleleAsUnit=T)), col="blue", space=0, 
          xlab="Allele counts",  
          main=paste(levels(pop(pop)),sum(table(glSum(pop, alleleAsUnit=T))),"SNPs", 
                     sep=" ")) 
})  
dev.off() 

par(mfrow=c(1,1)) 

###run glPcaFast script 
y<-aa.genlight
pca.1 <- glPcaFast(aa.genlight, nf=300) 
# proportion of explained variance by first three axes 
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis 
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis 
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis 
# save fig 
pdf ("PCA_all_SNPs_ax12.pdf", width=14, height=7) 
col <- funky(5) 
g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, 
              col=transp(col,.6),  
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, 
              pgrid.draw =F, plot = FALSE) 
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = 
                 list(box = list(draw = FALSE),  
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE) 
ADEgS(c(g1, g2), layout = c(1, 2)) 
dev.off() 

#playing with PCA
summary(pca.1)
scores<-data.frame(pca.1$scores)
summary(scores)
meta$PC1log<-scores$PC1
meta$PC2log<-scores$PC2
meta$PC3log<-scores$PC3

names(meta)

m1<-glm(PC2log~RUN+island+YEAR, data= meta)
m2<-glm(PC2log~island+YEAR, data= meta)
m3<-glm(PC2log~RUN+YEAR, data= meta)
m4<-glm(PC2log~RUN+island, data= meta)
m5<-glm(PC2log~RUN, data= meta)
m6<-glm(PC2log~island, data= meta)
m7<-glm(PC2log~YEAR, data=meta)
AIC(m1,m2,m3,m4,m5,m6,m7)
logLik(m1)
logLik(m2)
logLik(m3)
logLik(m4)
logLik(m5)
logLik(m6)
logLik(m7)

library(car)
Anova(m6)

### K-means clustering

grp <- find.clusters(aa.genlight, max.n.clust=40, parallel=FALSE)
2
grp <- find.clusters(aa.genlight, max.n.clust=3000, glPca = pca.1, perc.pca = 
                       100, n.iter=1e6, n.start=1000)

write.table(grp$grp, file="grouping_Kmeans_all.txt", sep="\t", quote=F, 
            col.names=F) 

### Calculate Nei's distances between individuals/pops 
aa.D.ind <- stamppNeisD(y, pop = FALSE)  # Nei's 1972 distance between indivs 
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst") # export matrix - for SplitsTree 
aa.D.pop <- stamppNeisD(y, pop = TRUE)   # Nei's 1972 distance between pops 
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance.phy.dst") # export matrix - for SplitsTree 

summary(aa.D.ind)
### Calculate pairwise Fst among populations 
aa.genlight@ploidy <- as.integer(ploidy(y)) 
aa.fst<-stamppFst(y, nboots = 1, percent =95, nclusters=3) 
#modify the matrix for opening in SplitsTree 
aa.fst.sym <- aa.fst 
aa.fst.sym[upper.tri(aa.fst.sym)] <- t(aa.fst.sym)[upper.tri(aa.fst.sym)]   
# add upper triangle 
aa.fst.sym[is.na(aa.fst.sym)] <- 0                                  
#replace NAs with zero 
stamppPhylip(aa.fst.sym, file="ALL_aa.pops_pairwise_Fst.phy.dst")   # export matrix - for SplitsTree

colnames(aa.D.ind) <- rownames(aa.D.ind)    
pdf(file="Neis_dist_heatmap.pdf", width=10, height=10) 
heatmap.2(aa.D.ind, trace="none", cexRow=0.4, cexCol=0.4) 
dev.off() 

# plot and save NJ tree 
plot(nj(aa.D.ind)) 
write.tree(nj(aa.D.ind),file="NJ.Neis.dist.tree.tre") 

aa.genlight2 <- aa.genlight 
setPop(aa.genlight2)<-~RUN
#pop(aa.genlight2)<-substr(indNames(aa.genlight2),5,9)  # define populations as the AAXXX codes 

aa.D.pop2 <- stamppNeisD(aa.genlight2, pop = TRUE)     # Nei's 1972 distance between pops 
stamppPhylip(aa.D.pop2, file="aa.pops2_Neis_distance.phy.dst") # export matrix - for SplitsTree 
# create the dist objects used in analyses below 
colnames(aa.D.ind) <- rownames(aa.D.ind)    
aa.D.ind.dist<-as.dist(aa.D.ind, diag=T) 
attr(aa.D.ind.dist, "Labels")<-rownames(aa.D.ind)          # name the rows of a matrix   
colnames(aa.D.pop2) <- rownames(aa.D.pop2)  
aa.D.pop.dist<-as.dist(aa.D.pop2, diag=T) 
attr(aa.D.pop.dist, "Labels")<-rownames(aa.D.pop2)          # name the rows of a matrix   

pops <- as.factor(pop(aa.genlight2))                        # define populations 
groups <- as.factor(pop(aa.genlight))     # define groups 
# one-level AMOVA 
(res <- pegas::amova(aa.D.ind.dist ~ pops))           # one-level AMOVA, default nperm=1000 
(res <- pegas::amova(aa.D.ind.dist ~ groups))
# hierarchical AMOVA  
(res <- pegas::amova(aa.D.ind.dist ~ pops/groups))    # hierarchical AMOVA

coords <- read.csv ("pop_coords.txt", sep ="\t")     # tab-separated file for all pops 
xy.coords.only<- subset(coords, select=c("lat","lon")) 
Dgeo <- dist(xy.coords.only) 

library(ggmap) 
map <- get_map(location =  c(lon = -53, lat = 50), zoom = 6) 
mapPoints <- ggmap(map) + geom_point(data = coords, aes(x = lon, y = lat, 
                                                        colour="blue")) + geom_text(data = coords, aes(x = lon, y = lat,label = pop, 
                                                                                                       colour = "red"), size = 4, vjust = 0, hjust = -0.5) 
mapPoints 

IBD <- mantel.randtest(Dgeo,aa.D.ind.dist) 
IBD 
plot(Dgeo,aa.D.ind.dist, pch=20,cex=.5) 
abline(lm(aa.D.ind.dist~Dgeo)) 
Dgeo

#plot and check for denser areas in the plot indicating sub-groups 
library(MASS) 
dens <- kde2d(Dgeo,aa.D.ind.dist, n=300, lims=c(-1, 16, 0, 0.08)) 
summary(dens)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red")) 
plot(Dgeo, aa.D.ind.dist, pch=20,cex=.5) 
image(dens, col=(myPal(300)), add=TRUE) 
abline(lm(aa.D.ind.dist~Dgeo)) 
title("Correlation of Genetic and Geographical distances") 
