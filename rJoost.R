###############################
#
# Housekeeping, load librarys and functions
#
###############################
library("phyloseq")
library("ggplot2")
library('rhdf5')
#library('biom')
library('vegan')
library("ggplot2")
library("gplots")
library("DESeq2")
#library("vsn")
library("ape")
library("vegan")
library("plyr")
library("dplyr")
library("reshape2")
library("igraph")
library("ggnetwork")
library("intergraph")
library("gridExtra")
library('structSSI')
library("ade4")
library("knitr")
#library("BiocStyle")
library("caret")

.cran_packages <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny",
                    "miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest",
                    "vegan", "plyr", "dplyr", "ggrepel", "nlme",
                    "reshape2","devtools", "PMA", "structSSI", "ade4",
                    "igraph", "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("phyloseq", "genefilter", "impute")

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}

.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

#tol color schemes
tol1qualitative=c("#4477AA")
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

generate_matrix <- function(x){
  indptr  = x$sample$matrix$indptr+1
  indices = x$sample$matrix$indices+1
  data    = x$sample$matrix$data
  nr = length(x$observation$ids)
  
  counts = sapply(2:length(indptr),function(i){
    x = rep(0,nr)
    seq = indptr[i-1]:(indptr[i]-1)
    x[indices[seq]] = data[seq]
    x
  })
  rownames(counts) = x$observation$ids
  colnames(counts) = x$sample$ids
  # I wish this next line wasn't necessary
  lapply(1:nrow(counts),function(i){
    counts[i,]
  })
}
generate_metadata <- function(x){
  metadata = x$metadata
  metadata = lapply(1:length(x$ids),function(i){
    id_metadata = lapply(metadata,function(j){
      if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
      else{ j[i] }
    })
    list(id = x$ids[i],metadata=id_metadata)
  })
  return(metadata)
}
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}
read_hdf5_biom<-function(file_input){
  x = h5read(file_input,"/",read.attributes = TRUE)
  data = generate_matrix(x)
  rows = generate_metadata(x$observation)
  columns = generate_metadata(x$sample)
  shape = c(length(data),length(data[[1]])) # dim(data)
  # Experimental -- need to actually load these from file
  id = attr(x,"id")
  vs = attr(x,"format-version")
  format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
  format_url = attr(x,"format-url")
  type = "OTU table"
  #type=attr(x,"type")
  generated_by = attr(x,"generated-by")
  date = attr(x,"creation-date")
  matrix_type = "dense"
  matrix_element_type = "int"
  
  namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
            rows,columns,shape,data)
}

#geomean function for vst
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



#alpha diversity curve function
set.seed(42)


calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity, check.names=FALSE), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity', check.names=FALSE)
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth' ,.progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

theme_set(theme_bw())

###Loading data
setwd("D:/Rfiles/16s")
y = read_hdf5_biom("bac_arch.biom")
tree = 'otus.tre'
library(biomformat)
dataRaw0 = import_biom(biom(y), tree)

sample_sums(dataRaw0) ## to count the number of reads?
dataRaw0 = prune_samples(sample_sums(dataRaw0)>=1000,dataRaw0) #remove below 1000
dataRaw0 <- subset_samples(dataRaw0, Dup == "N") #remove duplicates
sample_sums(dataRaw0) ##again count


#convert factor
sample_data(dataRaw0)$AgeStage <- factor(sample_data(dataRaw0)$AgeStage)
sample_data(dataRaw0)$Bb <- factor(sample_data(dataRaw0)$Bb)
sample_data(dataRaw0)$Site <- factor(sample_data(dataRaw0)$Site)
sample_data(dataRaw0)$Design <- factor(sample_data(dataRaw0)$Design)
sample_data(dataRaw0)$Design2 <- factor(sample_data(dataRaw0)$Design2)
sample_data(dataRaw0)$Run <- factor(sample_data(dataRaw0)$Run)
sample_data(dataRaw0)$Dup <- factor(sample_data(dataRaw0)$Dup)
sample_data(dataRaw0)$Source <- sample_data(dataRaw0)$Bb [sample_data(dataRaw0)$Bb]
sample_data(dataRaw0)$Source <- factor(sample_data(dataRaw0)$Source)
sample_data(dataRaw0)  
  
write.csv(otu_table(dataRaw0.glom.6), file= "b.csv")
write.csv(tax_table(dataRaw0.glom.6),file= "a.csv")

summary(sample_data(dataRaw0))


###############################
#
# Create gloms
#
###############################

dataRaw0.glom.6 = tax_glom(dataRaw0, 'Rank6', NArm = FALSE)
dataRaw0.glom.5 = tax_glom(dataRaw0.glom.6, 'Rank5', NArm = FALSE)
dataRaw0.glom.4 = tax_glom(dataRaw0.glom.5, 'Rank4', NArm = FALSE)
dataRaw0.glom.3 = tax_glom(dataRaw0.glom.4, 'Rank3', NArm = FALSE)
dataRaw0.glom.2 = tax_glom(dataRaw0.glom.3, 'Rank2', NArm = FALSE)




###############################
#
# Preliminary manual insepction of dataset
# Remove very low abundance stuff which are more then likely noise
# and add nothing of value to our analysis at this point in time 
#
###############################

rank_names(dataRaw0.glom.6)
table(tax_table(dataRaw0.glom.6)[, "Rank2"], exclude = NULL)

prevdf = apply(X = otu_table(dataRaw0.glom.6),
               MARGIN = ifelse(taxa_are_rows(dataRaw0.glom.6), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(dataRaw0.glom.6),
                    tax_table(dataRaw0.glom.6))

View(plyr::ddply(prevdf, "Rank2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence),mean(df1$TotalAbundance),sum(df1$TotalAbundance))}))

###############################
#
# Misc function for calculating the tabel in the manuscript
#
###############################

genfac = factor(rownames(tax_table(dataRaw0)))
#genfac = factor(tax_table(dataRaw0)[, "Rank2"])
gentab = apply(otu_table(dataRaw0), MARGIN = 2 , function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1
#gentabsum = apply(gentab < 2, 2, sum)
gentabsum = apply(gentab > 1, 2, sum)
View(as.data.frame(gentabsum))

reads<-sampleSums(dataRaw0)
View(as.data.frame(reads))

###############################
#
# Create a shadow dataset, but here we glommerate everything
# and remove the two soil samples
#
###############################

dataTick = dataRaw0

#optional denoise-ing here
dataTick = subset_samples(dataTick, Source == "Soil") #remove duplicates
sample_sums(dataTick)

dataTick.glom.6 = tax_glom(dataTick, 'Rank6', NArm = FALSE)
dataTick.glom.5 = tax_glom(dataTick.glom.6, 'Rank5', NArm = FALSE)
dataTick.glom.4 = tax_glom(dataTick.glom.5, 'Rank4', NArm = FALSE)
dataTick.glom.3 = tax_glom(dataTick.glom.4, 'Rank3', NArm = FALSE)
dataTick.glom.2 = tax_glom(dataTick.glom.3, 'Rank2', NArm = FALSE)

###############################
#
# Create a shadow dataset, but here we glommerate everything
# and remove the two soil samples
#
###############################

dataNoSoil = dataRaw0

#optional denoise-ing her
dataTick = subset_samples(dataTick, Source == "Soil") #remove duplicates
sample_sums(dataTick)

otu_table(dataNoSoil)
dataNoSoil =  prune_taxa(taxa_sums(dataTick) > 0, dataNoSoil) 
sample_sums(dataNoSoil)

dataTick.glom.6 = tax_glom(dataTick, 'Rank6', NArm = FALSE)
dataTick.glom.5 = tax_glom(dataTick.glom.6, 'Rank5', NArm = FALSE)
dataTick.glom.4 = tax_glom(dataTick.glom.5, 'Rank4', NArm = FALSE)
dataTick.glom.3 = tax_glom(dataTick.glom.4, 'Rank3', NArm = FALSE)
dataTick.glom.2 = tax_glom(dataTick.glom.3, 'Rank2', NArm = FALSE)

###############################
#
# Figure 1A: Alpha diversity curves
#
###############################

#alphaDenoise <- filter_taxa(dataRaw, function(x) sum(x > 1) > (1), TRUE)
alphaDenoise =  prune_taxa(taxa_sums(dataTick) > 1, dataTick) 

rarefaction_curve_data <- calculate_rarefaction_curves(alphaDenoise, c('Observed', 'Simpson'), rep(c(1, 10, 100, 1000, 2000, 5000, 10000), each = 10))

#fix R's retardness here
rarefaction_curve_data$Sample <- gsub("X", "", rarefaction_curve_data$Sample)
rarefaction_curve_data$Sample <- gsub("\\.", "-", rarefaction_curve_data$Sample)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(alphaDenoise)), by.x = 'Sample', by.y = 'row.names')

ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Site,
    group = SampleID
  )) + geom_line()  + geom_point(size=1) + facet_wrap(
    facets = ~ Measure,
    scales = 'free_y'
  )

###############################
#
# Figure 1B: Alpha diversity based on complete specimens
#
###############################

#we can probbably better use the more then once, in at least 2 samples
#it supressed the spurious OTU's a bit.. 

alphaDenoise <- filter_taxa(dataTick, function(x) sum(x > 2) > (1), TRUE)
#alphaDenoise =  prune_taxa(taxa_sums(dataRaw) > 1, dataRaw) 


plot_richness(alphaDenoise, measures=c("Observed", "Simpson", "invsimpson"), x="Bb", sortby="chao1") + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot_richness(alphaDenoise, measures=c("Observed", "Simpson", "invsimpson"), x="Site" , sortby="Chao1") + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


frame = as.data.frame(estimate_richness(alphaDenoise,  measures=c("observed", "simpson", "invsimpson")))
frame$Bb = sample_data(alphaDenoise)$Bb
frame$Site = sample_data(alphaDenoise)$Site
View(frame)
summary(frame)
lm<-lm(InvSimpson~Bb, data=frame)
anova(lm)

plot_richness(alphaDenoise, measures=c("Observed", "Simpson", "invsimpson"), x="Site",  sortby="Simpson") + geom_boxplot()
plot_richness(alphaDenoise, measures=c("Observed", "Simpson", "invsimpson"), x="SampleID", sortby="Simpson")

#in text
estimate_richness(alphaDenoise, measures=c("invsimpson"))

###############################
#
# Figure 2A: Oridnation NMDS on Jensen-Shannon divergence
# using log transformed data
#
# Figure 2B: Distances plotted on hierarchical clustering graph
#
###############################

denoise <- filter_taxa(dataTick, function(x) sum(x > 10) > (1), TRUE)


sample_data(denoise)$Site <- factor(sample_data(denoise)$Site)

dataRawLog <- transform_sample_counts(denoise, function(x) log(1 + x))
out.wuf.log <- ordinate(dataRawLog, method = "MDS", distance = "jsd")
evals <- out.wuf.log$values$Eigenvalues

plot_ordination(dataRawLog, out.wuf.log, color = "Bb",
                shape = "Site") +
  scale_shape_manual(values=16:nlevels(sample_data(dataRawLog)$Site)) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size=3) +
  labs(col = "Bb", shape = "Site")

distance<-phyloseq::distance(dataRawLog, method='jsd')
clustered<-hclust(distance, method="ward.D2")

plot(clustered, hang = -1)
plot(as.phylo(clustered), type = "unrooted")

summary(sample_data(dataRawLog))
sample_data(dataRawLog)$DesignSite<-ifelse(sample_data(dataRawLog)$Site=="GULL",16,17)
sample_data(dataRawLog)$DesignBb<-ifelse(sample_data(dataRawLog)$Bb=="Pos",16,17)
                                           
library(ggtree)
ggtree(as.phylo(clustered),branch.length ="height",layout = "fan") + 
geom_tippoint(col=sample_data(dataRawLog)$DesignBb, shape=sample_data(dataRawLog)$DesignSite, size=3)




##
library("DESeq2")

dataDESeq = subset_samples(dataRaw0.glom.6)
sample_sums(dataDESeq)

diagdds = phyloseq_to_deseq2(dataDESeq, ~ Site)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = T)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataDESeq)[rownames(sigtab), ], "matrix"))
head(sigtab)

summary(res)
dim(sigtab)

summary(sigtab)

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Rank2, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Rank2), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Rank6, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Rank6), levels=names(x))
ggplot(sigtab, aes(x=Rank6, y=log2FoldChange, color=Rank2)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  facet_grid(.~Rank2, scales = "free_x")

ggplot(sigtab, aes(y=Rank6, x=log2FoldChange, color=Rank2)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  facet_grid(Rank2~.,scales="free", space = "free") +
  theme(strip.text.y = element_text(angle = 0))


###############################
#
# Composition heatmap for sponge individuals
#
###############################


###############################
#
# Attempt to combine sponges in same figure 
#
###############################

#grab samplse and denoise 
dataIndDenoised <- filter_taxa(dataTick, function(x) sum(x > 1) > (1), TRUE)
dataIndDenoisedLog <- transform_sample_counts(dataIndDenoised, function(x) log(1 + x))

#create a dendrogram for next to the plot
dataIndDistance<-phyloseq::distance(dataIndDenoisedLog, method='jsd')
dataIndClustered<-hclust(dataIndDistance, method="ward.D2")

#phylum 
choL2 = dataTick.glom.2
choL2.norm <- transform_sample_counts(choL2, function(x) (x/sum(x)) * 100)
gp.cho <- prune_taxa(names(sort(taxa_sums(choL2.norm), decreasing=T)),choL2.norm)
topList.cho = prune_taxa(names(sort(taxa_sums(gp.cho), decreasing=T)[1:4]),gp.cho)
taxonomy.cho = as.data.frame(tax_table(topList.cho))
theme_set(theme_bw())
plot_heatmap(topList.cho, sample.label="Bb")
cnames.cho <- paste(taxonomy.cho$Rank2, sep =' ')
cnames.cho = factor(append(as.character(cnames.cho),"Other"))
rnames.cho <- paste(sample_data(choL2)$Species, '-', sample_data(choL2)$Label, sep="")
dat.cho = as.data.frame(otu_table(topList.cho))
tdat.cho = t(dat.cho)
#dat.cho$Other = 100 - rowSums(dat.cho)
tdat.cho
dat2.cho = as.matrix(dat.cho)

summary(dat2.cho)
rownames(dat2.cho)
cnames.cho
rownames(dat2.cho) <- cnames.cho
#colnames(dat2.cla) <- cnames.cla
dat3.cho <- dat2.cho[,order(colSums(dat2.cho), decreasing = TRUE)]
#rownames(dat2.cho) <- cnames.cho
dat3dataframe.cho = as.data.frame(dat3.cho)
col_idx.cho <- grep("Other", names(dat3dataframe.cho))

dat3.cho <- dat3dataframe.cho[, c(col_idx.cho, (1:ncol(dat3dataframe.cho))[-col_idx.cho])]
dat3.cho <- dat3dataframe.cho[, c((1:ncol(dat3dataframe.cho))[-col_idx.cho], col_idx.cho)]
dat3.cho = as.matrix(dat3.cho)

my_palette <- colorRampPalette(c("yellow", "red", "blue"))(n = 100)
heatmap(dat2.cho)
heatmap(dat2.cho)
heatmap.2(dat2.cho,
          main = " ", # heat map title
          symm = FALSE,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          #dendrogram="row",# only draw a row dendrogram
          #Rowv=as.dendrogram(dataIndClustered),
          #Colv="NA", 
          key.xlab="Sequence abundance (%)",   srtCol=45)  


