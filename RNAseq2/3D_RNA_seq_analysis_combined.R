#####################################################################################
"""Script for analysis of RNA-seq salmon output (transcript abundance data). Salmon was run on two RNAseq datasets with barley Introgression lines with differing microbiota compliments. This is a custom version of the 3D RNA-seq pipeline by Wenbin Guo (https://github.com/wyguo/ThreeDRNAseq)
Max Coulter made the following changes/additions:
- General changes to allow for simultaneous analysis of two RNA-seq experiments
- Use of UMAP as well as PCA for dimensionality reduction and data visualisation
- Changing edgeR linear model to include date as a block effect
- Retrieval of raw p values for use in volcano plots (but significance is set by adjusted p values). It is recommended to display raw p expressed_genes_both_adj_p_values
- Volcano plots
- Go enrichment analyses, GO slim analyses, bubble plots
- Jitterplots for displaying individual gene cpm

Authors: Max Coulter and Wenbin Guo"""

candidate_genes <- "/cluster/db/mecoulter/RNAseq2/RNAquant_analysis_data/genes_in_region_bart2.csv"

gene_annotation <- "/cluster/db/mecoulter/BaRT2v18/BaRT_2_18_annotation_genes.txt"






#install.packages(statmod)

library(statmod)
library(tximport)
library(edgeR)
library(limma)
library(RUVSeq)
library(eulerr)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
#devtools::install_github("jokergoo/ComplexHeatmap")
#library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(plyr)
#library(pca3d)
library(dplyr)
#install.packages("magrittr")
#install.packages("tibble")
#install.packages("tidyverse")
library(tidyverse)
library(magrittr)
library(tibble)
#install.packages("beeswarm")
library("beeswarm")
library(umap)
library(VennDiagram)
#install.packages("circlize")
library(circlize)

library(UpSetR)


options(stringsAsFactors = F)

source("/cluster/db/mecoulter/RNAseq2/3D_source_code.R")#3D functions, from Wenbin Guo


DDD.data <- list()
################################################################################
##----->> Set folders to read and save results
data.folder <- file.path(getwd(),'RNAquant_analysis_data') # for .RData format results
result.folder <- file.path(getwd(),'RNAquant_analysis_result') # for 3D analysis results in csv files
figure.folder <- file.path(getwd(),'RNAquant_analysis_figures')# for figures
report.folder <- file.path(getwd(),'RNAquant_analysis_report')

DDD.data$data.folder <- data.folder
DDD.data$result.folder <- result.folder
DDD.data$figure.folder <- figure.folder
DDD.data$report.folder <- report.folder

if(!file.exists(data.folder))
  dir.create(path = data.folder,recursive = T)
if(!file.exists(result.folder))
  dir.create(path = result.folder,recursive = T)
if(!file.exists(figure.folder))
  dir.create(path = figure.folder,recursive = T)
if(!file.exists(report.folder))
  dir.create(path = report.folder,recursive = T)

### Set the input data folder
##----->> folder of input files
input.folder <- 'salmoncombined'
quant.folder <- 'quants'

quant_method <- 'salmon' # abundance generator
tximport_method <- 'lengthScaledTPM' # method to generate expression in tximport

################################################################################
##----->> parameters for data pre-processing
### has sequencign replicates?
has_srep <- T

### parameter for low expression filters
cpm_cut <- 2
cpm_samples_n <- 3

### data normalisation parameter
norm_method <- 'TMM' ## norm_method is one of 'TMM','RLE' and 'upperquartile'

################################################################################
##----->> parameters for 3D analysis
pval_adj_method <- 'BH'
pval_cut <- 0.01
l2fc_cut <- 1

deltaPS_cut <- 0.1
DAS_pval_method <- 'F-test'



################################################################################
##----->> Meta table includes sample information, e.g. conditions, bio-reps, seq-reps, abundance paths, etc.
#metatable <- read.csv(file.path(getwd(),'metadata-rhizosphere_RNA-Seq.csv'))
metatable <- read.csv(file.path(getwd(),'metadata-rhizosphere_RNA-Seqexperiment2only.csv'))
##select the columns of experimental design
factor_col <-'Tissue'
brep_col <- 'brep'
srep_col <- 'srep'
quant_col <- 'quant.files'
experiment <- 'experiment'
sample_date <- 'sampling.date'

##arrange information in the metatable
metatable$label <- as.vector(interaction(metatable[,factor_col]))
metatable$sample.name <- as.vector(interaction(metatable[,c(factor_col,brep_col,srep_col,experiment)]))
metatable$quant.folder <- file.path(quant.folder,metatable$quant.files, ifelse(quant_method=='salmon','quant.sf','abundance.h5'))

##----->> Transcript-gene association mapping
mapping <-read.csv(file.path(getwd(),'transcripts_genes_BaRT2v18.csv'))

mapping <- data.frame(as.matrix(mapping),stringsAsFactors = F)
rownames(mapping) <- mapping$TXNAME


################################################################################
##----->> Generate gene expression
##
setwd(input.folder)
txi_genes <- tximport(metatable$quant.folder,dropInfReps = T,
	type = quant_method, tx2gene = mapping,
	countsFromAbundance = tximport_method)

## give columnnames to the datasets
colnames(txi_genes$counts) <-
colnames(txi_genes$abundance) <-
colnames(txi_genes$length) <-metatable$sample.name

## save the data
write.csv(txi_genes$counts,file=paste0(result.folder,'/counts_genes.csv'))
write.csv(txi_genes$abundance,file=paste0(result.folder,'/TPM_genes.csv'))
save(txi_genes,file=paste0(data.folder,'/txi_genes.RData'))

################################################################################
##----->> Generate transcripts expression
txi_trans<- tximport(metatable$quant.folder,
	type = quant_method, tx2gene = NULL,
	countsFromAbundance = tximport_method,
	txOut = T,dropInfReps = T)

## give colunames to the datasets
colnames(txi_trans$counts) <-
  colnames(txi_trans$abundance) <-
  colnames(txi_trans$length) <-metatable$sample.name

## save the data
write.csv(txi_trans$counts,file=paste0(result.folder,'/counts_trans.csv'))
write.csv(txi_trans$abundance,file=paste0(result.folder,'/TPM_trans.csv'))
save(txi_trans,file=paste0(data.folder,'/txi_trans.RData'))

################################################################################
##extract gene and transcript read counts
genes_counts <- txi_genes$counts
trans_counts <- txi_trans$counts
trans_TPM <- txi_trans$abundance


##If no sequencing replicates, genes_counts and trans_counts remain the same by

################################################################################

if(has_srep){
  idx <- paste0(metatable$label,'.',metatable[,brep_col],'.',metatable$experiment)
  genes_counts <- sumarrays(genes_counts,group = idx)
  trans_counts <- sumarrays(trans_counts,group = idx)
  metatable_new <- metatable[metatable[,srep_col]==metatable[,srep_col][1],]
} else {
  metatable_new <- metatable
}

metatable_new$count_name <- colnames(trans_counts)





##Now we do neeed the 3D RNA-seq package...

#install.packages("remotes")
#remotes::install_github("wyguo/ThreeDRNAseq")




##----->> Do the filters

target_high <- low.expression.filter(abundance = trans_counts,
	mapping = mapping,
	abundance.cut = cpm_cut,
	sample.n = cpm_samples_n,
	unit = 'counts',
	Log=F)



##save expressed genes and transcripts
save(target_high,file=paste0(data.folder,'/target_high.RData'))

################################################################################
##----->> Mean-variance plot
## transcript level

counts.raw = trans_counts[rowSums(trans_counts>0)>0,]
counts.filtered = trans_counts[target_high$trans_high,]
mv.trans <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.trans$fit.raw
fit.filtered <- mv.trans$fit.filtered
mv.trans.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (transcript level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (transcript level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.trans.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Transcript mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.trans.plot()
dev.off()

################################################################################
## gene level
counts.raw = genes_counts[rowSums(genes_counts>0)>0,]
counts.filtered = genes_counts[target_high$genes_high,]
mv.genes <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.genes$fit.raw
fit.filtered <- mv.genes$fit.filtered
mv.genes.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (gene level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (gene level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.genes.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Gene mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.genes.plot()
dev.off()

#pdf(file = paste0(figure.folder,'/Gene mean-variance trend.pdf'),
    #width = 25/2.54,height = 12/2.54)
#mv.genes.plot()
#dev.off()
####PCA plot########
################################################################################
##----->> trans level


data2pca <- trans_counts[target_high$trans_high,]
dge <- DGEList(counts=data2pca)
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
groups <- metatable_new[,brep_col] ## colour on biological replicates

g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

#g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA Bio-reps.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: genotypes',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

png(filename = paste0(figure.folder,'/Transcript PCA genotypes.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()
#pdf(file = paste0(figure.folder,'/Transcript PCA Bio-reps.pdf'),
    #width = 15/2.54,height = 13/2.54)
#print(g)
#dev.off()

#Colour on experiment number

groups <- metatable_new$experiment ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: experiment',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

png(filename = paste0(figure.folder,'/Transcript PCA experiemnt.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

#Now see effect of sampling date (extra for root RNA-Seq)
groups <- metatable_new$sampling.date ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: date',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

png(filename = paste0(figure.folder,'/Transcript PCA date.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

#Now look just at experiemnt 2
second_experiment_metatable_new <- subset(metatable_new,metatable_new$experiment==2,droplevels=TRUE)

se_trans_counts <- trans_counts[,second_experiment_metatable_new$count_name]

data2pca <- se_trans_counts[target_high$trans_high,]
dge <- DGEList(counts=data2pca)
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')


groups <- second_experiment_metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: genotypes',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

png(filename = paste0(figure.folder,'/Transcript PCA genotypes second experiment.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()
#####

pca <-prcomp(metabo[,-1],scale.=TRUE)
#gr  <-factor(metabo[,1])summary(gr)
#pca3d(pca,group=gr)
#snapshotPCA3d(file="first_plot.png")

####
##################################################
##--average expression plot
groups <- metatable_new[,brep_col]

data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'Transcript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

#pdf(file = paste0(figure.folder,'/Transcript PCA average expression.pdf'),
    #width = 15/2.54,height = 13/2.54)
#print(g)
#dev.off()


################################################################################
##----->> genes level
data2pca <- genes_counts[target_high$genes_high,]
dge <- DGEList(counts=data2pca)
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots

groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Bio-reps.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: genotypes',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA genotypes.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

#Experiment

groups <- metatable_new$experiment ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: experiment',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA experiment.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()

#NOw try with sample date

groups <- metatable_new$sampling.date ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: date',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = T)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA date.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(g)
dev.off()


#pdf(file = paste0(figure.folder,'/Gene PCA Bio-reps.pdf'),
    #width = 15/2.54,height = 13/2.54)
#print(g)
#dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

#pdf(file = paste0(figure.folder,'/Gene PCA average expression.pdf'),
    #width = 15/2.54,height = 13/2.54)
#print(g)
#dev.off()
#######################
##----->> genes level
data2pca <- genes_counts[target_high$genes_high,]
dge <- DGEList(counts=data2pca)
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--treatment plots

#groups <- metatable_new[,brep] ## colour on treatment
groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()


#################################UMAP###############################


source("https://github.com/crj32/M3C/blob/master/R/umap.R")#This doesn't work. So code from repo pasted below...

umap <- function(mydata, labels=FALSE, printres=FALSE, seed=FALSE, axistextsize = 18,
                 legendtextsize = 18, dotsize = 5, textlabelsize = 4, legendtitle = 'Group',
                 controlscale = FALSE, scale = 1, low = 'grey', high = 'red', 
                 colvec = c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen", 
                            "violetred", "orange", "midnightblue", "grey31", "black"),
                 printheight = 20, printwidth = 22, text = FALSE){
  
  ## basic error handling
  
  if ( controlscale == TRUE && class(labels) %in% c( "character", "factor") && scale %in% c(1,2) ) {
    stop("when categorical labels, use scale=3")
  }
  if ( controlscale == TRUE && class(labels) %in% c( "numeric") && scale %in% c(3) ) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if ( controlscale == FALSE && scale %in% c(2,3) ) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c('character','factor')){
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- 'Unknown'
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c('character','factor')){
    warning("there is NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- 'Unknown'
  }
  
  ##
  
  message('***UMAP wrapper function***')
  message('running...')
  
  if (seed != FALSE){
    set.seed(seed)
  }
  
  if (labels[1] == FALSE && text[1] == FALSE){
    
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(colour='skyblue', size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize))+
      scale_colour_manual(values = colvec)
    
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAP.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] != FALSE && text[1] == FALSE){ ##### KEY
    
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    
    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + #scale_colour_distiller(palette = "Spectral")
          scale_colour_gradient(low=low, high=high)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) +
          scale_colour_manual(values = colvec)
      }
    }else{
      p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) + 
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)
    }
    
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAPlabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] != FALSE && text[1] != FALSE){ ##### KEY
    
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    scores$label <- text
    
    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + #scale_colour_distiller(palette = "Spectral")
          scale_colour_gradient(low=low, high=high)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) +
          scale_colour_manual(values = colvec)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }
    }else{
      p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) + 
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    }
    
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAPlabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] == FALSE && text[1] != FALSE){
    
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    scores$label <- text
    
    p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + 
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = colvec) + geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAP.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }
  
  message('done.')
  
  return(p)
}

metatable$genotype.experiment <- as.vector(interaction(metatable[,c(factor_col,experiment)]))
metatable_new$genotype.experiment <- as.vector(interaction(metatable_new[,c(factor_col,experiment)]))
metatable_new$experiment <- as.factor(metatable_new$experiment)

#UMAP at gene level:
data2pca <- genes_counts[target_high$genes_high,]

#Experiment 1


first_experiment_metatable_new <- subset(metatable_new,metatable_new$experiment==1,droplevels=TRUE)
fe_gene_counts <- data2pca[,first_experiment_metatable_new$count_name]

UMAP <- umap(fe_gene_counts,labels=first_experiment_metatable_new$label,colvec=c("green","blue","purple"))

png(filename = paste0(figure.folder,'/Gene UMAP first experiment genotype2.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()

#Dates

UMAP <- umap(fe_gene_counts,labels=first_experiment_metatable_new$sampling.date,colvec=c("black","dark grey","light grey"))

png(filename = paste0(figure.folder,'/Gene UMAP first experiment date.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()



#Now look just at experiemnt 2
second_experiment_metatable_new <- subset(metatable_new,metatable_new$experiment==2,droplevels=TRUE)

se_gene_counts <- data2pca[,second_experiment_metatable_new$count_name]



UMAP <- umap(se_gene_counts,labels=second_experiment_metatable_new$label,colvec=c("blue","green","purple"))

png(filename = paste0(figure.folder,'/Gene UMAP second experiment genotype2.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()


png(filename = paste0(figure.folder,'/Gene UMAP second experiment genotype2.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(umap(se_gene_counts,labels=second_experiment_metatable_new$label,colvec=c("blue","green","purple")))
dev.off()



#Date

UMAP <- umap(se_gene_counts,labels=second_experiment_metatable_new$sampling.date,colvec=c("black","dark grey","light grey"))

png(filename = paste0(figure.folder,'/Gene UMAP second experiment date.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()


#Now with both experiments:



UMAP <- umap(data2pca,labels=metatable_new$genotype.experiment,colvec=c("green","yellow","red","blue","dark blue","dark green"))

png(filename = paste0(figure.folder,'/Gene UMAP both experiments genotype.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()

#Date
UMAP <- umap(data2pca,labels=metatable_new$sampling.date,colvec=c("black","dark grey","light grey"))

png(filename = paste0(figure.folder,'/Gene UMAP both experiments date.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()

#Experiment number

UMAP <- umap(data2pca,labels=metatable_new$experiment,colvec=c("purple","green"))

png(filename = paste0(figure.folder,'/Gene UMAP both experiments experiment.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()

#Microbial phenotype

UMAP <- umap(data2pca,labels=metatable_new$phenotype,colvec=c("blue","orange"))

png(filename = paste0(figure.folder,'/Gene UMAP both experiments phenotype.png'),
    width = 20/2.54,height = 17/2.54,units = 'in',res = 300)
print(UMAP)
dev.off()




####Normalisation#################
################################################################################
##----->> trans level
dge <- DGEList(counts=trans_counts[target_high$trans_high,],
               group = metatable_new$label,
               genes = mapping[target_high$trans_high,])
trans_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(trans_dge,file=paste0(data.folder,'/trans_dge.RData'))

################################################################################
##----->> genes level
dge <- DGEList(counts=genes_counts[target_high$genes_high,],
               group = metatable_new$label)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(genes_dge,file=paste0(data.folder,'/genes_dge.RData'))

################################################################################
##----->> distribution plot
sample.name <- paste0(metatable_new$label,'.',metatable_new[,brep_col])
condition <- metatable_new$label

###--- trans level
data.before <- trans_counts[target_high$trans_high,]
data.after <- counts2CPM(obj = trans_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Transcript expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

#pdf(file = paste0(figure.folder,'/Transcript expression distribution.pdf'),
    #width = 20/2.54,height = 20/2.54)
#do.call(grid.arrange,g)
#dev.off()

###--- genes level
data.before <- genes_counts[target_high$genes_high,]
data.after <- counts2CPM(obj = genes_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Gene expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

#pdf(file = paste0(figure.folder,'/Gene expression distribution.pdf'),
    #width = 20/2.54,height = 20/2.54)
#do.call(grid.arrange,g)
#dev.off()
###############Data information##############################################
RNAseq_info <- data.frame(
  Description=c('Raw transcripts',
                'Raw genes',
                'Samples',
                'Samples after merging seq-reps',
                'Condition of interest',
                'CPM cut-off',
                'Min samples to CPM cut-off',
                'Expressed transcripts',
                'Expressed genes'),
  Number=c(length(mapping$TXNAME),
           length(unique(mapping$GENEID)),
           nrow(metatable),
           nrow(metatable_new),
           length(unique(metatable$label)),
           cpm_cut,
           cpm_samples_n,
           length(target_high$trans_high),
           length(target_high$genes_high))
)
DDD.data$RNAseq_info <- RNAseq_info

RNAseq_info

#####Make contrast groups#####################
################################################################################
##----->> pair-wise contrast groups
contrast <- c('HID19.1-Barke.1','HID56.1-Barke.1','HEB_124_17.2-Barke.2','HEB_124_52.2-Barke.2','HEB_124_52.2-HEB_124_17.2')



batch.effect <- NULL ## if has no batch effects




#Oringinal

design <- condition2design(condition = metatable_new$genotype.experiment,
                           batch.effect = batch.effect)


################################################################################

##----->> edgeR glmQL pipeline
genes_3D_stat_control <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glmQL',
                                 adjust.method = pval_adj_method)




## edgeR
DE_pipeline = 'glmQL'



design <- model.matrix(~0 + genotype.experiment + sampling.date, data=metatable_new)

idx <- grep('genotype.experiment',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('genotype.experiment')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design

idx <- grep('sampling.date',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('sampling.date')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design


colnames(design) <- make.names(colnames(design))


##----->> edgeR glmQL pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glmQL',
                                 adjust.method = pval_adj_method)



## save results
DDD.data$genes_3D_stat <- genes_3D_stat
DDD.data$genes_3D_stat_control <- genes_3D_stat_control


################################################################################
##----->> generate deltaPS
deltaPS <- transAbundance2PS(transAbundance =txi_trans$abundance[target_high$trans_high,],
                             PS = NULL,
                             contrast = contrast,
                             condition = metatable$genotype.experiment,
                             mapping = mapping[target_high$trans_high,])

DDD.data$PS <- PS <- deltaPS$PS
DDD.data$deltaPS <- deltaPS <- deltaPS$deltaPS


################################################################################
##----->> DAS genes,DE and DTU transcripts
#batch.effect <- genes_batch$W
#batch.effect <- NULL ## if has no batch effects
#design <- condition2design(condition = metatable_new$label,
                           #batch.effect = batch.effect)

################################################################################



trans_3D_stat <- edgeR.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 method = 'glmQL',
                                 adjust.method = pval_adj_method)


## save results
DDD.data$trans_3D_stat <- trans_3D_stat

################################################################################
##----->> Summary DE genes
DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval_cut,
                                         log2FC=l2fc_cut))

DE_genes_control <- summaryDEtarget(stat = genes_3D_stat_control$DE.stat,
                              cutoff = c(adj.pval=pval_cut,
                                         log2FC=l2fc_cut))

DDD.data$DE_genes <- DE_genes

################################################################################
## summary DAS genes, DE and DTU trans
##----->> DE trans
DE_trans <- summaryDEtarget(stat = trans_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval_cut,
                                         log2FC=l2fc_cut))
DDD.data$DE_trans <- DE_trans

##----->> DAS genes
if(DAS_pval_method=='F-test') {
  DAS.stat <- trans_3D_stat$DAS.F.stat
} else {
  DAS.stat <- trans_3D_stat$DAS.Simes.stat
}

lfc <- genes_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DAS_genes <- summaryDAStarget(stat = DAS.stat,
                                lfc = lfc,
                                cutoff=c(pval_cut,deltaPS_cut))
DDD.data$DAS_genes <- DAS_genes

##----->> DTU trans
lfc <- trans_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DTU_trans <- summaryDAStarget(stat = trans_3D_stat$DTU.stat,
                                lfc = lfc,cutoff = c(adj.pval=pval_cut,
                                                     deltaPS=deltaPS_cut))
DDD.data$DTU_trans <- DTU_trans

################################################################################
## save csv
write.csv(DE_genes,file=paste0(result.folder,'/DE genes.csv'),row.names = F)
write.csv(DAS_genes,file=paste0(result.folder,'/DAS genes.csv'),row.names = F)
write.csv(DE_trans,file=paste0(result.folder,'/DE transcripts.csv'),row.names = F)
write.csv(DTU_trans,file=paste0(result.folder,'/DTU transcripts.csv'),row.names = F)


################################################################################
##----->> target numbers
DDD_numbers <- summary3Dnumber(DE_genes = DE_genes,
                                 DAS_genes = DAS_genes,
                                 DE_trans = DE_trans,
                                 DTU_trans=DTU_trans,
                                 contrast = contrast)
DDD_numbers

#control
DDD_numbers_control <- summary3Dnumber(DE_genes = DE_genes_control,
                                 DAS_genes = DAS_genes,
                                 DE_trans = DE_trans,
                                 DTU_trans=DTU_trans,
                                 contrast = contrast)
DDD_numbers_control
#The model (that includes date as a batch effect) has the biggest effect on HID19-Barke, with many more DE genes




write.csv(DDD_numbers,file=paste0(result.folder,'/DE DAS DTU numbers.csv'),
          row.names = F)
DDD.data$DDD_numbers <- DDD_numbers
################################################################################
##----->> DE vs DAS
DEvsDAS_results <- DEvsDAS(DE_genes = DE_genes,
                           DAS_genes = DAS_genes,
                           contrast = contrast)
DEvsDAS_results
DDD.data$DEvsDAS_results <- DEvsDAS_results
write.csv(DEvsDAS_results,file=paste0(result.folder,'/DE vs DAS gene number.csv'),
          row.names = F)


################################################################################
##----->> DE vs DTU
DEvsDTU_results <- DEvsDTU(DE_trans = DE_trans,
                           DTU_trans = DTU_trans,
                           contrast = contrast)
DEvsDTU_results
DDD.data$DEvsDTU_results <- DEvsDTU_results
write.csv(DEvsDTU_results,file=paste0(result.folder,'/DE vs DTU transcript number.csv'),row.names = F)


#Make plots
################################################################################
##----->> DE genes
idx <- factor(DE_genes$contrast,levels = contrast)
targets <-  split(DE_genes,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down-regulated','up-regulated'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DE genes',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DE genes up and down regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

#pdf(paste0(figure.folder,'/DE genes up and down regulation numbers.pdf'),
    #width = length(contrast)*5/2.54,10/2.54)
#print(g.updown)
#dev.off()

################################################################################
##----->> DE trans
idx <- factor(DE_trans$contrast,levels = contrast)
targets <-  split(DE_trans,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down-regulated','up-regulated'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DE trans',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DE transcripts up and down regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

#pdf(paste0(figure.folder,'/DE transcripts up and down regulation numbers.pdf'),
    #width = length(contrast)*5/2.54,10/2.54)
#print(g.updown)
#dev.off()
################################################################################
##----->> DE genes
#Eurler diagram (Author Max Coulter)


contrast2plot <- c("HID19-Barke","HID56-Barke")










HID19 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID19-Barke.1",drop=TRUE)))
HID56 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID56-Barke.1",drop=TRUE)))




#palette using black
cbbPalette <- c("#F0E442",'#D55E00','#0072B2')



venn.diagram(
  x = list(HID19$target,HID56$target),
  category.names = c("Int19 - Barke" , "Int56 - Barke"),
  filename = 'DE_genes_venn.png',
  output=TRUE,
  # Output features
        imagetype="svg" ,
        height = 1000 ,
        width = 1000 ,
        resolution = 300,
        compression = "lzw",
          # Circles
        lwd = 0.5,
        fill = cbbPalette,
        cat.cex = 1,
        cat.pos = c(-27, 27),
        cat.dist = c(0.035, 0.035)

)
dev.off()





#3 way venn diagram for experiment 2:

#Try first with the control....

HEB52_HEB17 <- subset(DE_genes_control,DE_genes_control$contrast=="HEB_124_52.2-HEB_124_17.2",drop=TRUE)
HEB52_Barke <- subset(DE_genes_control,DE_genes_control$contrast=="HEB_124_52.2-Barke.2",drop=TRUE)
HEB17_Barke <- subset(DE_genes_control,DE_genes_control$contrast=="HEB_124_17.2-Barke.2",drop=TRUE)

#Intersects:
HEB52_HEB17_HEB52_Barke <-intersect(HEB52_HEB17$target,HEB52_Barke$target)#51

HEB52_Barke_HEB17_Barke <- intersect(HEB52_Barke$target, HEB17_Barke$target)#55
HEB52_HEB17_HEB17_Barke <- intersect(HEB52_HEB17$target, HEB17_Barke$target)#55


all_intersect <- intersect(intersect(HEB52_HEB17, HEB52_Barke), HEB17_Barke) #0

length(setdiff(setdiff(HEB52_Barke$target,HEB52_Barke_HEB17_Barke), HEB52_HEB17_HEB52_Barke))

length(setdiff(setdiff(HEB17_Barke$target,HEB52_Barke_HEB17_Barke), HEB52_HEB17_HEB17_Barke))

length(setdiff(setdiff(HEB52_HEB17$target,HEB52_HEB17_HEB17_Barke), HEB52_HEB17_HEB52_Barke))






cbbPalette <- c('#E69F00','#D55E00','#0072B2')
cbbPalette <- c("#F0E442","#56B4E9",'purple')
venn.diagram(
  x = list(HEB52_HEB17$target,HEB52_Barke$target,HEB17_Barke$target),
  category.names = c("HEB_124_52 - Barke", "HEB_124_17 - Barke", "HEB_124_52 - HEB_124_17"),
  filename = paste0(figure.folder, '/DE_genes_venn_experiment2_control.svg'),
  output=TRUE,
  # Output features
        imagetype="svg" ,
        height = 100 ,
        width = 100 ,
        resolution = 300,
        compression = "lzw",
          # Circles
        lwd = 0.5,
        alpha = c(0.3,0.3,0.3),
        cat.cex = 0.7,
        cat.pos = c(-15, 15, 180),
        cat.dist = c(0.035, 0.035, 0.035)
        #fill = cbbPalette

)
dev.off()


venn.diagram(
  x = list(HEB52_HEB17$target,HEB52_Barke$target),
  category.names = c("HEB_124_52 - Barke", "HEB_124_52 - HEB_124_17"),
  filename = paste0(figure.folder, '/DE_genes_venn_experiment2_control2.svg'),
  output=TRUE,
  # Output features
        imagetype="svg" ,
        height = 100 ,
        width = 100 ,
        resolution = 300,
        compression = "lzw",
          # Circles
        lwd = 0.5,
        alpha = c(0.3,0.3),
        cat.cex = 0.7,
        cat.pos = c(-15, 15),
        cat.dist = c(0.035, 0.035)
        #fill = cbbPalette

)
dev.off()




HEB52_HEB17 <- subset(DE_genes,DE_genes$contrast=="HEB_124_52.2-HEB_124_17.2",drop=TRUE)
HEB52_Barke <- subset(DE_genes,DE_genes$contrast=="HEB_124_52.2-Barke.2",drop=TRUE)
HEB17_Barke <- subset(DE_genes,DE_genes$contrast=="HEB_124_17.2-Barke.2",drop=TRUE)

#Intersects:
length(intersect(HEB52_HEB17$target,HEB52_Barke$target))#51

length(intersect(HEB52_Barke$target, HEB17_Barke$target))#55

HEB17_Barke_HEB52_HEB17 <- intersect(HEB17_Barke$target, HEB52_HEB17$target)

all_intersect <- intersect(intersect(HEB52_HEB17, HEB52_Barke), HEB17_Barke) #0





cbbPalette <- c('#E69F00','#D55E00','#0072B2')
cbbPalette <- c("#F0E442","#56B4E9",'purple')
venn.diagram(
  x = list(HEB52_HEB17$target,HEB52_Barke$target,HEB17_Barke$target),
  category.names = c("HEB_124_52 - Barke", "HEB_124_17 - Barke", "HEB_124_52 - HEB_124_17"),
  filename = paste0(figure.folder, '/DE_genes_venn_experiment2.svg'),
  output=TRUE,
  # Output features
        imagetype="svg" ,
        height = 100 ,
        width = 100 ,
        resolution = 300,
        compression = "lzw",
          # Circles
        lwd = 0.5,
        fill = cbbPalette,
        alpha = c(0.3,0.3,0.3),
        cat.cex = 0.7,
        cat.pos = c(-15, 15, 180),
        cat.dist = c(0.035, 0.035, 0.035)

)
dev.off()


#####UpSet plots#####
png(filename= paste0(figure.folder,"/upset_figure.png"), width = 900, height = 900, units = "px") # or other device
x = list("124_52 - 124_17" = HEB52_HEB17$target, "124_52 - Barke" = HEB52_Barke$target, "124_17 - Barke" = HEB17_Barke$target)
print({upset(fromList(x), order.by = "freq", point.size = 3.5, line.size = 3,
mainbar.y.label = "", sets.x.label = "No. genes", text.scale = c(3.3, 3.3, 3, 3, 3, 3))})
dev.off()






plotVolcano <- function(raw_p,
						lgfc,
						title='',
						col0='black',col1='blue',
						col2='vermillion',
						size=1.5,
						ylab='-log10(p)',
						xlab='log2FC',
						pval_adj_method=pval_adj_method,
						candidate_gene_path=candidate_genes,
            gg_force=8){
	HID19_raw<-as.data.frame(cbind(raw_p,lgfc))
	adjusted<-p.adjust(HID19_raw$raw_p, method = pval_adj_method)
	HID19_raw<-cbind(HID19_raw,adjusted)
	significance <- NULL
	for(i in 1:nrow(HID19_raw)) {
		if (HID19_raw$adjusted[i] <= 0.01 & (HID19_raw$lgfc[i] >= 1 | HID19_raw$lgfc[i] <= -1 )) {
			if (HID19_raw$lgfc[i] >= 1) {
			significance = c(significance,"significant up")
			} else {
			significance = c(significance,"significant down")
			}
		} else {
			significance = c(significance,"not significant")
		}
	}
	HID19_raw<-cbind(HID19_raw,significance)
	candidate_genes <- as.data.frame(read.csv(candidate_gene_path,header = TRUE))
	HID19_candidates <- subset(HID19_raw,row.names(HID19_raw) %in% candidate_genes$gene,drop=TRUE)
	#Get gene functions
	print(paste("Number of expressed candidate genes in interval ",length(row.names(HID19_candidates)),sep = ""))
	gene_functions <- NULL
	for(i in 1:nrow(candidate_genes)) {
		if (candidate_genes$gene[i] %in% row.names(HID19_candidates)) {
			print(candidate_genes$gene[i])
			gene_functions = c(gene_functions,candidate_genes$function.[i])
		}
	}
	HID19_candidates_functions <- cbind(HID19_candidates,gene_functions)
	write.csv(HID19_candidates_functions,file=paste0(result.folder,'/HID19_candidates_functions.csv'),row.names = T)
	png(filename = paste("volcano_",title,"_rawp_ggplot.png",sep = ""),width = 900, height = 900, units = "px", pointsize = 12,bg = "white")
	print({
	g <- ggplot(data = HID19_raw,aes(x=lgfc,y=-log10(raw_p))) + geom_point(aes(colour=significance),size=size) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_x_continuous(breaks = pretty(HID19_raw$lgfc, n = 10)) + scale_y_continuous(breaks = pretty(-log10(HID19_raw$raw_p), n = 10)) + geom_label_repel(data=HID19_candidates_functions,aes(x=lgfc,y=-log10(raw_p),label=gene_functions,color = significance),size = 4, show.legend=FALSE, fontface = 'bold',box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'purple',arrow = arrow(length = unit(0.01, 'npc')),force=8,nudge_x=4,nudge_y=0.5) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
	})
	dev.off()
	#Subset candidate list to only include those that are significant
	HID19_candidates_functions_sig <- subset(HID19_candidates_functions,(HID19_candidates_functions$significance=='significant up'|HID19_candidates_functions$significance=='significant down'))
	png(filename = paste("volcano_",title,"_rawp_signamesonly_ggplot.png",sep = ""),width = 9000, height = 9000, units = "px", pointsize = 12,bg = "white",res=1000)
	print({
		g <- ggplot(data = HID19_raw,aes(x=lgfc,y=-log10(raw_p))) + geom_point(aes(colour=significance),size=size) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_x_continuous(breaks = pretty(HID19_raw$lgfc, n = 10)) + scale_y_continuous(breaks = pretty(-log10(HID19_raw$raw_p), n = 10)) + geom_label_repel(data=HID19_candidates_functions_sig,aes(x=lgfc,y=-log10(raw_p),label=gene_functions,color = significance),size = 5, show.legend=FALSE, fontface = 'bold',box.padding = unit(0.35, "lines"),point.padding = unit(0.2, "lines"),segment.color = 'black',arrow = arrow(length = unit(0.01, 'npc')),force=gg_force,nudge_x=6,nudge_y=3) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
	})
	dev.off()
}




#output expressed genes from both RNAseq experiments

write.csv(genes_3D_stat$DE.pval,file=paste0(result.folder,'/expressed_genes_both_adj_p_values.csv'))


#HID19
raw_p<-genes_3D_stat$DE.pval[,1]
lgfc<-genes_3D_stat$DE.lfc[,1]


plotVolcano(raw_p,lgfc,title='HID19-Barke',col0='black',col1='#E69F00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HID56
raw_p<-genes_3D_stat$DE.pval[,2]
lgfc<-genes_3D_stat$DE.lfc[,2]

plotVolcano(raw_p,lgfc,title='HID56-Barke',col0='black',col1='#D55E00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HID19
raw_p<-genes_3D_stat$DE.pval[,1]
lgfc<-genes_3D_stat$DE.lfc[,1]

#Above are not raw p values... try again

raw_p<-genes_3D_stat$DE.rawpval[,1]
lgfc<-genes_3D_stat$DE.lfc[,1]


plotVolcano(raw_p,lgfc,title='Int19-Barke_corrected',col0='black',col1='#E69F00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HID56
raw_p<-genes_3D_stat$DE.rawpval[,2]
lgfc<-genes_3D_stat$DE.lfc[,2]

plotVolcano(raw_p,lgfc,title='Int56-Barke_corrected',col0='black',col1='#D55E00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HEB_124_17.2-Barke.2
raw_p<-genes_3D_stat$DE.rawpval[,3]
lgfc<-genes_3D_stat$DE.lfc[,3]

plotVolcano(raw_p,lgfc,title='HEB_124_17-Barke_corrected',col0='black',col1='#D55E00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HEB_124_52.2-Barke.2
raw_p<-genes_3D_stat$DE.rawpval[,4]
lgfc<-genes_3D_stat$DE.lfc[,4]

plotVolcano(raw_p,lgfc,title='HEB_124_52-Barke_corrected',col0='black',col1='#D55E00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method)

#HEB_124_52.2-HEB_124_17.2
raw_p<-genes_3D_stat$DE.rawpval[,5]
lgfc<-genes_3D_stat$DE.lfc[,5]

plotVolcano(raw_p,lgfc,title='HEB_124_52-HEB_124_17',col0='black',col1='#E69F00',col2='#0072B2',size=1.5,ylab='-log10(p)',xlab='log2FC',pval_adj_method=pval_adj_method,gg_force=15)

#Plot logfc of genes according to position
annotation <- read.csv("/cluster/db/mecoulter/BaRT2v18/BaRT_2_18_annotation_genes.txt",sep = "\t", header = TRUE)
rownames(annotation) <- annotation$BaRTv2.gene




plotPosition <- function(adj_p,
						lgfc,
						title='',
						col0='black',col1='blue',
						col2='vermillion',
						size=1.5,
            chromosome,
						annotation){
lgfc2 <- lgfc
lgfc <-as.data.frame(cbind(adj_p,lgfc))
lgfc$gene_name <- rownames(lgfc)
#gene2term <- gene2term %>% filter(gene_id %in% genes_all)
lgfc_merged <- merge(lgfc, annotation, by.x = "gene_name", by.y = "BaRTv2.gene", all.x = TRUE, all.y = FALSE)

significance <- NULL
	for(i in 1:length(adj_p)) {
		if (adj_p[i] <= 0.01 & (lgfc2[i] >= 1 | lgfc2[i] <= -1 )) {
			if (lgfc2[i] >= 1) {
			significance = c(significance,"significant up")
			} else {
			significance = c(significance,"significant down")
			}
		} else {
			significance = c(significance,"not significant")
		}
	}

lgfc_merged$significance <- significance

lgfc_merged_3H <- subset(lgfc_merged, lgfc_merged$Chromosome==chromosome,drop=TRUE)
#NOw plot according to start position

ylab<-'Log2FC'
xlab<-paste0('Position on ', chromosome,' (Mbp)')
b_i<-100000000#break interval
png(filename = paste(chromosome,title,"_lgfc_start.png",sep = ""),width = 9000, height = 4000, units = "px", pointsize = 12,bg = "white",res=1000)
	print({
		g <- ggplot(data = lgfc_merged_3H,aes(x=Start,y=lgfc)) + geom_point(aes(colour=significance),size=size) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_y_continuous(breaks=pretty(lgfc_merged_3H$lgfc, n = 10)) + scale_x_continuous(breaks=c(0,1*b_i,2*b_i,3*b_i,4*b_i,5*b_i,6*b_i), labels=c(0,100,200,300,400,500,600)) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
	})
dev.off()
}

#HEB_124_52.2-Barke.2
adj_p<-genes_3D_stat$DE.pval[,4]
lgfc<-genes_3D_stat$DE.lfc[,4]


plotPosition(adj_p,lgfc,title='HEB_124_52.2-Barke.2',col0='black',col1='#E69F00',col2='#0072B2',chromosome="chr3H",annotation=annotation)

#HEB124_52 - HEB_124_17
adj_p<-genes_3D_stat$DE.pval[,5]
lgfc<-genes_3D_stat$DE.lfc[,5]

plotPosition(adj_p,lgfc,title='HEB_124_52.2-HEB_124_17.2',col0='black',col1='#E69F00',col2='#56B4E9',chromosome="chr3H",annotation=annotation)

#All chromosomes
for (i in 1:7){
  chromosome_name <- paste0("chr",i,"H")
 plotPosition(adj_p,lgfc,title='HEB_124_52.2-HEB_124_17.2',col0='black',col1='#E69F00',col2='#56B4E9',chromosome=chromosome_name,annotation=annotation)
}




#NOw HEB124_52 - HEB_124_17 subset in Hid124_52-Barke comparison
col0='black'
col1='#E69F00'
col2='#56B4E9'
adj_p52<-genes_3D_stat$DE.pval[,4]
lgfc52<-genes_3D_stat$DE.lfc[,4]
chromosome = "chr3H"
title='HEB_124_52.2-HEB_124_17.2_subset'

merged_52 <- as.data.frame(cbind(adj_p52,lgfc52,adj_p,lgfc))

genessig_52_17 <- subset(DE_genes,DE_genes$contrast=="HEB_124_52.2-HEB_124_17.2",drop=TRUE)
genessig_52_Barke <- subset(DE_genes,DE_genes$contrast=="HEB_124_52.2-Barke.2",drop=TRUE)


significance <- NULL
	for(i in 1:nrow(merged_52)) {
		if ((adj_p[i] <= 0.01 & (lgfc[i] >= 1 | lgfc[i] <= -1 )) & (adj_p52[i] <= 0.01 & (lgfc52[i] >= 1 | lgfc52[i] <= -1 ))) {
			if (lgfc[i] >= 1 & lgfc52[i] >= 1 ){
			significance = c(significance,"significant up")
			} else {
        if (lgfc[i] <= -1 & lgfc52[i] <= -1 ){
			    significance = c(significance,"significant down")
        } else{
          significance = c(significance,"not significant")
          print("A wierd one...")
        }
			}
		} else {
			significance = c(significance,"not significant")
		}
	}

merged_52$significance <- significance

merged_52$gene_name <- rownames(merged_52)
#gene2term <- gene2term %>% filter(gene_id %in% genes_all)
merged_522 <- merge(merged_52, annotation, by.x = "gene_name", by.y = "BaRTv2.gene", all.x = TRUE, all.y = FALSE)


for (i in 1:7){
  chromosome <- paste0("chr",i,"H")
  lgfc_merged_3H <- subset(merged_522, merged_522$Chromosome==chromosome,drop=TRUE)
#NOw plot according to start position
 ylab='Log2FC'
 xlab=paste0('Position on ', chromosome, " (Mbp)")
 b_i<-100000000#break interval
 png(filename = paste0(figure.folder,"/",chromosome,title,"_lgfc_start.png"),width = 9000, height = 4000, units = "px", pointsize = 12,bg = "white", res = 1000)
	print({
		g <- ggplot(data = lgfc_merged_3H,aes(x=Start,y=lgfc)) + geom_point(aes(colour=significance),size=size) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_y_continuous(breaks = pretty(lgfc_merged_3H$lgfc, n = 10)) + scale_x_continuous(breaks=c(0,1*b_i,2*b_i,3*b_i,4*b_i,5*b_i,6*b_i), labels=c(0,100,200,300,400,500,600)) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
	})
 dev.off()
}



#Just experiment 2 results, no experiment 1 (not taking date into model)

#NOw HEB124_52 - HEB_124_17 subset in Hid124_52-Barke comparison
col0='black'
col1='#E69F00'
col2='#56B4E9'
adj_p52<-genes_3D_stat_control$DE.pval[,5]
lgfc52<-genes_3D_stat_control$DE.lfc[,5]
chromosome = "chr3H"
title='HEB_124_52.2-HEB_124_17.2_subset_2only'

merged_52 <- as.data.frame(cbind(adj_p52,lgfc52))

#Get significant genes of interest
HEB52_HEB17 <- DE_genes_control %>% filter(contrast == "HEB_124_52.2-HEB_124_17.2")
HEB52_Barke <- DE_genes_control %>% filter(contrast == "HEB_124_52.2-Barke.2")

HEB52_HEB17_HEB52_Barke <-intersect(HEB52_HEB17$target,HEB52_Barke$target)#51


#Add significance column with information whether up or down regulated

significance <- list()
all_genes <-rownames(merged_52)
for(i in 1:nrow(merged_52)) {
  if (all_genes[i] %in% HEB52_HEB17_HEB52_Barke){
    if (lgfc52[i] >= 1 ){
      significance[length(significance) + 1] <- "significant up"
    }else {
      if (lgfc52[i] <= -1){
       significance[length(significance) + 1] <- "significant down"
      } else {
        print("problem!")
      }
    }
  } else{
    significance[length(significance) + 1] <- "not significant"
  }
}





merged_52$significance <- unlist(significance)

merged_52$gene_name <- rownames(merged_52)
#gene2term <- gene2term %>% filter(gene_id %in% genes_all)
merged_522 <- merge(merged_52, annotation, by.x = "gene_name", by.y = "BaRTv2.gene", all.x = TRUE, all.y = FALSE)


for (i in 1:7){
  chromosome <- paste0("chr",i,"H")
  lgfc_merged_3H <- subset(merged_522, merged_522$Chromosome==chromosome,drop=TRUE)
#NOw plot according to start position
 ylab='Log2FC'
 xlab=paste0('Position on ', chromosome, " (Mbp)")
 b_i<-100000000#break interval
 png(filename = paste0(figure.folder,"/",chromosome,title,"_lgfc_start.png"),width = 9000, height = 4000, units = "px", pointsize = 12,bg = "white", res = 1000)
	print({
		g <- ggplot(data = lgfc_merged_3H,aes(x=Start,y=lgfc52)) + geom_point(aes(colour=significance),size=size) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_y_continuous(breaks = pretty(lgfc_merged_3H$lgfc52, n = 10)) + scale_x_continuous(breaks=c(0,1*b_i,2*b_i,3*b_i,4*b_i,5*b_i,6*b_i), labels=c(0,100,200,300,400,500,600)) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
	})
 dev.off()
}







merged19_52_sig <- subset(merged19_52, (merged19_52$significance == "significant up") | (merged19_52$significance == "significant down"), drop=TRUE)



plotPosition(adj_p,lgfc,title='HEB_124_52.2-HEB_124_17.2_subset',col0='black',col1='#E69F00',col2='#56B4E9',chromosome="chr3H",annotation=annotation)

####Plot logfc of both comparisons:

contrast2plot <- c("HID19-Barke","HID56-Barke")

HID19<-genes_3D_stat$DE.lfc[,1]
HID56<-genes_3D_stat$DE.lfc[,2]

logfc <- as.data.frame(cbind(HID19,HID56))

HID19 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID19-Barke",drop=TRUE)))
HID56 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID56-Barke",drop=TRUE)))

both <- as.data.frame(as.vector(subset(HID19,HID19$target %in% HID56$target,drop=TRUE)))

`%notin%` <- Negate(`%in%`)
onlyHID19 <- as.data.frame(as.vector(subset(HID19,HID19$target %notin% HID56$target,drop=TRUE)))
onlyHID56 <- as.data.frame(as.vector(subset(HID56,HID56$target %notin% HID19$target,drop=TRUE)))

significance <- NULL
for(i in 1:nrow(logfc)) {
		if (row.names(logfc)[i] %in% both$target) {
			significance = c(significance,"significant both")
		} else {
      if (row.names(logfc)[i] %in% onlyHID19$target) {
			  significance = c(significance,"significant Int19 only")
      } else {
        if (row.names(logfc)[i] %in% onlyHID56$target) {
			    significance = c(significance,"significant Int56 only")
        } else {
          significance = c(significance,"not significant")
        }
      }
    }
}


logfc <- as.data.frame(cbind(logfc,significance))
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(filename = "Logfcs.png",width = 600, height = 600, units = "px", pointsize = 12,bg = "white",quality = 100)
	print({
g <- ggplot(data = logfc,aes(x=HID19,y=HID56)) + geom_point(aes(colour=significance),size=1) + theme_bw() + labs(x="log2fc Int19 - Barke",y="log2fc Int56 - Barke",title = "") + scale_color_manual(values=c('not significant'="#999999",'significant Int19 only'="#D55E00", 'significant Int56 only'="#E69F00",'significant both'="#009E73")) + scale_x_continuous() + scale_y_continuous(breaks = c(-15,-10,-5,0,5,10))
})
dev.off()

##Ok now add labels of genes in interval

candidate_genes <- read.csv(candidate_gene_path,header = TRUE)
logfc_candidates <- subset(logfc,row.names(logfc) %in% candidate_genes$gene,drop=TRUE)
#Get gene functions
print(paste("Number of expressed candidate genes in interval ",length(row.names(logfc_candidates)),sep = ""))
gene_functions <- NULL
for(i in 1:nrow(candidate_genes)) {
	if (candidate_genes$gene[i] %in% row.names(logfc_candidates)) {
		print(candidate_genes$gene[i])
		gene_functions = c(gene_functions,candidate_genes$function.[i])
	}
}
logfc_candidates_functions <- cbind(logfc_candidates,gene_functions)

	logfc_candidates_functions_sig <- subset(logfc_candidates_functions,(logfc_candidates_functions$significance=='significant both'|logfc_candidates_functions$significance=='significant Int19 only'|logfc_candidates_functions$significance=='significant Int56 only'))

#FOr now just those genes with significance in both
logfc_candidates_functions_sig <- subset(logfc_candidates_functions,(logfc_candidates_functions$significance=='significant both'))

png(filename = "Logfcs_sigbothlabels.png",width = 600, height = 600, units = "px", pointsize = 12,bg = "white",quality = 100)
	print({
g <- ggplot(data = logfc,aes(x=HID19,y=HID56)) + geom_point(aes(colour=significance),size=1) + theme_bw() + labs(x="log2fc Int19 - Barke",y="log2fc Int56 - Barke",title = "") +  scale_x_continuous() + scale_y_continuous(breaks = c(-15,-10,-5,0,5,10)) + geom_label_repel(data=logfc_candidates_functions_sig,aes(x=HID19,y=HID56,label=gene_functions,color = significance),min.segment.length = 0, size = 3, show.legend=FALSE, fontface = 'bold',box.padding = unit(0.35, "lines"),point.padding = unit(0, "lines"),segment.color = 'black',arrow = arrow(length = unit(0.01, 'npc')),force=8,nudge_x=4,nudge_y=3) + scale_color_manual(values=c('not significant'="#999999",'significant Int19 only'="#999999", 'significant Int56 only'="#999999",'significant both'="#009E73"))
})
dev.off()


png(filename = "Logfcs.png",width = 600, height = 600, units = "px", pointsize = 12,bg = "white",quality = 100)
	print({
g <- ggplot(data = logfc,aes(x=HID19,y=HID56)) + geom_point(aes(colour=significance),size=1) + theme_bw() + labs(x="log2fc Int19 - Barke",y="log2fc Int56 - Barke",title = "") + scale_color_manual(values=c('not significant'="#999999",'significant Int19 only'="#999999", 'significant Int56 only'="#999999",'significant both'="#009E73")) + scale_x_continuous() + scale_y_continuous(breaks = c(-15,-10,-5,0,5,10))
})
dev.off()

# Pull out cpm plot of BART2_0-3H_12950,NBS-LRR like



#gene_name <- "BaRT2v18chr3HG123500"
gene_name <- "BaRT2v18chr3HG123500"


test2 <- t(t(genes_dge$counts[gene_name,]))
names <- rownames(test2)
colnames(test2) <- "CPM"
genotype_names <- NULL
genotypes <- NULL
experiment <- NULL
for (i in 1:length(names)) {
  split <- as.list(strsplit(names[i], "[.]")[[1]])
  name <- paste0(split[[1]], " (", split[[3]],")")
  genotype_name <- split[[1]]
  experiment <- c(experiment, split[[3]])
  genotype_names <- c(genotype_names,name)
  genotypes <- c(genotypes, genotype_name)
}


for_plot <- data.frame(CPM=test2,genotype=genotype_names, experiment_number=experiment, genotype_only= genotypes)

for_plot_2 <- subset(for_plot, for_plot$experiment=="2", drop=TRUE)



for_plot_2$genotype <- factor(for_plot_2$genotype , levels=c("Barke (2)", " HEB_124_17 (2)", " HEB_124_52 (2)"))

#colors <- c("#0072B2","#009E73","#CC79A7","#D55E00", "#56B4E9", "#E69F00")


colors <- c("#56B4E9","#E69F00","#0072B2")

png(filename = paste0(figure.folder, "/", gene_name, "_jitter.png"),width = 600, height = 600, units = "px", pointsize = 12,bg = "white")
	print({
g <- ggplot(for_plot_2, aes(x=genotype_only, y=CPM, colour=genotype_only)) +
    geom_boxplot(show.legend=FALSE) + theme_bw() + geom_jitter(width = 0.25, size = 5, aes(colour=genotype_only), show.legend=FALSE) + labs(x="Genotype",y="Normalised counts (CPM)",title = gene_name) + theme(axis.text.x=element_text(colour=colors,size=20),axis.text.y=element_text(size=25),
        axis.title=element_text(size=25),plot.title=element_text(size=25,hjust = 0.5)) + scale_color_manual(values=colors)
})
dev.off()

#####Clustering####

#First pick out the genes we need

#First get contrast groups
contrast2plot <- c("HID19-Barke","HID56-Barke")

HID19 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID19-Barke",drop=TRUE)))
HID56 <- as.data.frame(as.vector(subset(DE_genes,DE_genes$contrast=="HID56-Barke",drop=TRUE)))

all <- rbind(HID19,HID56)
#Up
all_up <- as.data.frame(as.vector(subset(all,all$up.down=="up-regulated",drop=TRUE)))
#Down
all_down <- as.data.frame(as.vector(subset(all,all$up.down=="down-regulated",drop=TRUE)))

all_up_5 <- as.data.frame(as.vector(subset(all_up,all_up$log2FC>=5,drop=TRUE)))

all_down_5 <- as.data.frame(as.vector(subset(all_down,all_down$log2FC<=-5,drop=TRUE)))
all_down_10 <- as.data.frame(as.vector(subset(all_down,all_down$log2FC<=-10,drop=TRUE)))

all_sig_genes <- unique(all$target)
#Now get counts for significant genes only
counts <- as.data.frame(genes_dge$counts)
counts_logcpm <- cpm(genes_dge, log=TRUE)
sig_genes_dge <-subset(counts, row.names(counts) %in% all_sig_genes, drop=TRUE)


#Heatmap:
#Now incorporate Wenbins code into what I have already done. Use my significant genes but use Wenbin's clustering - he uses raw counts instead of normalised counts
targets <- all_sig_genes

data2heatmap <- txi_genes$abundance[targets,]#Wenbin uses TPMs

data2heatmap <- sig_genes_dge
names <- colnames(data2heatmap)

genotype_names <- NULL
for (i in 1:length(names)) {
  split <- as.list(strsplit(names[i], "[.]")[[1]])
  genotype_names <- c(genotype_names,split[[1]])
}

data2plot <- rowmean(x = t(data2heatmap),
                     group = genotype_names,
                     reorder = F)

data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot)
hc <- hclust(hc.dist,method = "complete")
clusters <- cutree(hc, k = 8)


column_title <- paste0(length(targets),' DE genes')





#library(purr)

clusters <- reorderClusters(clusters = clusters,dat = hclust_matrix)

### save the target list in each cluster to result folder
x <- split(names(clusters),clusters)
x <- lapply(names(x),function(i){
  data.frame(Clusters=i,Targets=x[[i]])
})
x <- do.call(rbind,x)
colnames(x) <- c('Clusters','Targets')

###############################


#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cold="#009E73"
#hot=#E69F00"

col_fun = colorRamp2(c(-2, 0, 2), c("#009E73", "white", "#E69F00"))
#ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))


g <- Heatmap(as.matrix(data2plot), name = 'Z-scores',
             cluster_rows = TRUE,
             clustering_method_rows="complete",
             row_dend_reorder = T,
             show_row_names = FALSE,
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap_DE_genes.png'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54,
    units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()



#GO enrichment -

functionalEnrichment <- function(genes_all, genes_sel,  gene2term, term_info, gene2name = NULL,
                                 min_count=3, sig_limit=0.05) {
  #function author Marek Gerlinski, Univesity of Dundee
  #gene2term <- term_data$gene2term
  #term_info <- term_data$terms

  # select only terms represented in our gene set
  gene2term <- gene2term %>% filter(gene_id %in% genes_all)

  # all terms present in the selection
  terms <- gene2term %>%
    filter(gene_id %in% genes_sel) %>%
    pull(term_id) %>%
    unique()

  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)

  # empty line for missing terms
  na_term <- term_info %>% slice(1) %>% mutate_all(~NA)

  res <- map_dfr(terms, function(term) {
    info <- term_info %>% filter(term_id == term)
    # returns NAs if no term found
    if(nrow(info) == 0) info <- na_term %>% mutate(term_id = term)

    # all genes with the term
    tgenes <- gene2term %>% filter(term_id == term) %>% pull(gene_id)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)

    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)

    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow=2)
    ft <- fisher.test(fish, alternative = "greater")
    p <- as.numeric(ft$p.value)

    if(!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] %>% unname()

    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse=","),
        P = p
      )
    )
  }) %>%
    mutate(P = p.adjust(P, method="BH")) %>%
    filter(sel >= min_count & P <= sig_limit) %>%
    arrange(desc(enrich)) %>%
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))

  res
}


#Now import custom GO annotations
# Enrichment function
#
# Performs hypergeometric test.
# Input:
#   genes_sel and genes_all are vectors with gene IDs.
#   term_data: a list with two elements:
#      gene2term - data frame linking genes and terms, two columns: gene_id and term_id
#      terms - term description, data frame with at least two columns: term_id and decription; other columns will be kept
#   gene2name: named vector to change genes ids into gene names

setwd("/cluster/db/mecoulter/CEMartinez_RNAseq-209915706/results/GO/")

gene2term <- read.csv("BaRTv2.18_genes_go_ids.txt",sep="\t",header = FALSE)

terms <- read.csv("go_ids_descriptions.txt",sep="\t",header = TRUE)

gene2term <- gene2term %>% rename(gene_id = V1, term_id = V2)


term_data <- c(gene2term,terms)

genes_all <- as.vector(unique(genes_3D_stat$DE.stat$target))

#Remind oneself of the contrast groups
contrast_groups <- genes_3D_stat$DE.stat %>% pull(contrast) %>% unique()

DE_experiment1 <- genes_3D_stat$DE.stat %>% filter(contrast == "HID19.1-Barke.1" | contrast == "HID56.1-Barke.1")

sigall <- DE_experiment1 %>% filter(adj.pval <= 0.01) %>% filter(log2FC >= 1 | log2FC <= -1) %>% pull(target) %>% unique()

#DE_experiment1 %>% pull(contrast) %>% unique()

#Now pull out significant genes from both contrasts
Int19 <- DE_experiment1 %>% filter(contrast == "HID19.1-Barke.1") %>% filter(adj.pval <= 0.01) %>% filter(log2FC >= 1 | log2FC <= -1) %>% pull(target)

#Sanity check
(DDD_numbers %>% filter(contrast == "HID19.1-Barke.1") %>% pull("DE genes")) == length(Int19)



Int56 <- DE_experiment1 %>% filter(contrast == "HID56.1-Barke.1") %>% filter(adj.pval <= 0.01) %>% filter(log2FC >= 1 | log2FC <= -1) %>% pull(target)

#Sanity check
(DDD_numbers %>% filter(contrast == "HID56.1-Barke.1") %>% pull("DE genes")) == length(Int56)

sigboth <- DE_experiment1 %>% filter(target %in% Int19) %>% filter(target %in% Int56) %>% pull(target) %>% unique()



enriched_GOs <- functionalEnrichment(genes_all, genes_sel=sigboth, gene2term, terms)

write.csv(enriched_GOs, paste0(result.folder, "/", "19_56_sig_both_GO.csv"))

#This time will all 19_56 genes:

enriched_GOs <- functionalEnrichment(genes_all, genes_sel=sigall, gene2term, terms)

write.csv(enriched_GOs, paste0(result.folder, "/", "19_56_sig_all_GO.csv"))

#Actually just use GO slims:
gene2term <- read.csv("BaRTv2.18_gene_slim.txt",sep="\t",header = TRUE)

terms <- read.csv("goslim_ids_descriptions.txt",sep="\t",header = TRUE)


enriched_GOs <- functionalEnrichment(genes_all, genes_sel=sigboth, gene2term, terms)

write.csv(enriched_GOs, paste0(result.folder, "/", "19_56_sig_both_GO_slims.csv"))

enriched_GOs <- functionalEnrichment(genes_all, genes_sel=sigall, gene2term, terms)

write.csv(enriched_GOs, paste0(result.folder, "/", "19_56_sig_all_GO_slims.csv"))


genes_up <- as.vector(unique(all_up$target))#1245 genes with logfc >=1 in either group
genes_down <- as.vector(unique(all_down$target)) #740 genes with logfc <=-1 in either group

genes_up_5 <- as.vector(unique(all_up_5$target)) #56 genes with logfc >=5 in either group
genes_down_5 <- as.vector(unique(all_down_5$target))#160 genes with logfc <=-5 in either group
genes_down_10 <- as.vector(unique(all_down_10$target))
genes_down_8 <- as.vector(unique(all_down_8$target))
#First try with all significant genes
genes_sel <- all_sig_genes#All significant genes

enriched_GOs <- functionalEnrichment(genes_all, genes_sel, gene2term, terms)

enriched_GOs_up <- functionalEnrichment(genes_all, genes_up, gene2term, terms)
enriched_GOs_down <- functionalEnrichment(genes_all, genes_down, gene2term, terms)

enriched_GOs_up_5 <- functionalEnrichment(genes_all, genes_up_5, gene2term, terms)#Just ADP binding
enriched_GOs_down_5 <- functionalEnrichment(genes_all, genes_down_5, gene2term, terms)
enriched_GOs_down_8 <- functionalEnrichment(genes_all, genes_down_8, gene2term, terms)
enriched_GOs_down_10 <- functionalEnrichment(genes_all, genes_down_10, gene2term, terms)


#Get annotation...
annotation <- read.csv("BaRT_2_10_rRNAremoved_annotation_genes.txt",sep="\t")
annotation$target<-annotation$BaRT.2.gene

#Join together
annotation_up <- join(all_up,annotation,by="target",type="left")
annotation_down <- join(all_down,annotation,by="target",type="left")


#Write out the genes...

write.csv(annotation_up,"DEgenes_up.csv")
write.csv(annotation_down,"DEgenes_down.csv")





write.csv(enriched_GOs,"enriched_GOs_all_DEgenes.csv")
write.csv(enriched_GOs_up,"enriched_GOs_DEgenesup.csv")
write.csv(enriched_GOs_down,"enriched_GOs_DEgenesdown.csv")

write.csv(enriched_GOs_up_5,"enriched_GOs_DEgenesup5.csv")
write.csv(enriched_GOs_down_5,"enriched_GOs_DEgenesdown5.csv")

colnames(x) <- c('Clusters','Targets')
write.csv(x,"8_clusters.csv")
cluster_numbers <- unique(x$Clusters)

all_enriched <- data.frame()
for(i in 1:max(cluster_numbers)) {
  print(paste0("Working on cluster ",i))
  genes_s <-subset(x,x$Clusters == i,droplevels=TRUE)
  genes_sel <- genes_s$Targets
  enriched_GOs <- functionalEnrichment(genes_all, genes_sel, gene2term, terms)
  cluster <- rep(i,nrow(enriched_GOs))
  all_enriched <- rbind(all_enriched,cbind(enriched_GOs,cluster))
}


write.csv(all_enriched,"enriched_GOs_clusters.csv")


#Now import GO slims - created using owltools
gene2term <- read.csv("BaRTv2.18_gene_slim.txt",sep="\t",header = TRUE)

terms <- read.csv("goslim_ids_descriptions.txt",sep="\t",header = TRUE)

genes_all <- as.vector(unique(genes_3D_stat$DE.stat$target))

genes_up <- as.vector(unique(all_up$target))#1245 genes with logfc >=1 in either group
genes_down <- as.vector(unique(all_down$target)) #740 genes with logfc <=-1 in either group

genes_up_5 <- as.vector(unique(all_up_5$target)) #56 genes with logfc >=5 in either group
genes_down_5 <- as.vector(unique(all_down_5$target))#160 genes with logfc <=-5 in either group
genes_down_10 <- as.vector(unique(all_down_10$target))
genes_down_8 <- as.vector(unique(all_down_8$target))
#First try with all significant genes
genes_sel <- all_sig_genes#All significant genes


enriched_GOs_up <- functionalEnrichment(genes_all, genes_up, gene2term, terms)
enriched_GOs_down <- functionalEnrichment(genes_all, genes_down, gene2term, terms)

enriched_GOs_up_5 <- functionalEnrichment(genes_all, genes_up_5, gene2term, terms)



write.csv(enriched_GOs_up,"enriched_slims_DEgenesup.csv")
write.csv(enriched_GOs_down,"enriched_slims_DEgenesdown.csv")
write.csv(enriched_GOs_up_5,"enriched_slims_DEgenesup5.csv")

#Create bubble plots of GO enrichment

png(filename = "DE_slims_up_bubble.png",width = 600, height = 600, units = "px", pointsize = 12,bg = "white",quality = 100)
	print({
g <- ggplot(enriched_GOs_up, aes(x=sel, y=-log10(P), size=enrich)) + theme_bw() + geom_point(alpha=0.7,colour="#0072B2") + scale_size(range = c(2, 15),name="Enriched") + geom_label_repel(data=enriched_GOs_up,aes(x=sel,y=-log10(P),label=description),min.segment.length = 0, size = 3, show.legend=FALSE, fontface = 'bold',segment.color = 'black',force=8,nudge_x=2,nudge_y=2) + labs(x="Number of selected genes with GO slim",y="-log10(p)",title = "Significantly Enriched GO slims") + theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),plot.title=element_text(size=20,hjust = 0.5))
})
dev.off()


+ labs(x="Genotype",y="Normalised counts (CPM)",title = "BART2_0-3H_12950") + theme(axis.text.x=element_text(colour=c("#0072B2",'#E69F00',"#D55E00"),size=25),axis.text.y=element_text(size=25),
        axis.title=element_text(size=25),plot.title=element_text(size=25,hjust = 0.5)) + scale_color_manual(values=c("#0072B2",'#E69F00',"#D55E00"))



#GO slim analysis of 52 17 comparison

gene2term <- read.csv("/cluster/db/mecoulter/CEMartinez_RNAseq-209915706/results/GO/BaRTv2.18_gene_slim.txt",sep="\t",header = TRUE)

terms <- read.csv("/cluster/db/mecoulter/CEMartinez_RNAseq-209915706/results/GO/goslim_ids_descriptions.txt",sep="\t",header = TRUE)

genes_all <- as.vector(unique(genes_3D_stat$DE.stat$target))

genessig_52_17 <- subset(DE_genes,DE_genes$contrast=="HEB_124_52.2-HEB_124_17.2",drop=TRUE)


enriched_GOslim_52_17 <- functionalEnrichment(genes_all, genessig_52_17$target, gene2term, terms) #Nothing

#GO enrichment analysis:

gene2term <- as.data.frame(read.csv("BaRTv2.18_genes_go_ids.txt",sep="\t",header = FALSE))
gene2term$gene_id <- gene2term$V1
gene2term$term_id <- gene2term$V2

geneterms <- as.data.frame(cbind(gene2term$gene_id,gene2term$term_id))

terms <- as.data.frame(read.csv("go_ids_descriptions.txt",sep="\t",header = TRUE))



enriched_GOs_52_17 <- functionalEnrichment(genes_all, genessig_52_17$target, gene2term, terms) #Nothing

#Candidate gene information combined (52_17 comparison):

annotation <- as.data.frame(read.csv(gene_annotation,header = TRUE))

candidates <- as.data.frame(read.csv(candidate_genes,header = TRUE))

adj_p <- genes_3D_stat$DE.pval[,5]
raw_p <- genes_3D_stat$DE.rawpval[,5]
lgfc <- genes_3D_stat$DE.lfc[,5]


all_52_17 <- cbind(raw_p,adj_p,lgfc)
all_52_17 <- cbind(rownames(all_52_17), data.frame(all_52_17, row.names=NULL))
colnames(all_52_17)[1] <- "genes"


all_52_17_annotation <- merge(all_52_17, annotation, by.x = "genes", by.y = "BaRTv2.gene", all.x = TRUE, all.y = FALSE)
all_52_17_annotation_3H <- merge(candidates, all_52_17_annotation, by.x = "gene", by.y = "genes", all.x = FALSE, all.y = FALSE)

#Get rid of superfluous columns
keeps <- c("gene","Chromosome", "Start", "End","Strand" ,"raw_p", "adj_p", "lgfc", "function.","GO.IDs" )
all_52_17_annotation_3H_clean = all_52_17_annotation_3H[keeps]

write.csv(all_52_17_annotation_3H_clean,paste0(result.folder,"/candidate_gene_complete_info.csv"))




#Output list of the 50 DE genes that are of interest globally. Include information on position, RNA-seq, function (use dplyr syntax)



#BaRTv2 annotation file

BaRTv2 <- read_tsv(gene_annotation) %>% select(`BaRTv2 gene`, Chromosome, Start, End, Strand, `Pannzer annotation`) %>% rename(target = `BaRTv2 gene`)

locus3H_genes <- read_tsv(gene_annotation) %>% filter(Chromosome == "chr3H") %>% filter(Start >= 33181340) %>% filter(End <= 36970860) %>% select(`BaRTv2 gene`, Chromosome, Start, End, `Pannzer annotation`, `Coding potentiality`)

read_tsv(gene_annotation) %>% filter(Start > End) %>% nrow() == 0 # Test, should be True

write.csv(locus3H_genes, paste0(result.folder, "/", "all_3H_genes.csv"))



#Subset of 50 that are present in 124_52 - 124_17, 124_52 - Barke and absent from 124_17 - Barke
DE_genest <- as_tibble(DE_genes_control) #Convert to tibble so all things work ok

possible_contrasts <- DE_genest %>% select(contrast) %>% unique()

experiment2 <- DE_genest %>% filter(str_split_fixed(contrast,fixed("."), n = 3)[,3] == 2)

c52_17 <- experiment2 %>% filter(contrast == "HEB_124_52.2-HEB_124_17.2")

c52_Barke <- experiment2 %>% filter(contrast == "HEB_124_52.2-Barke.2")

c17_Barke <- experiment2 %>% filter(contrast == "HEB_124_17.2-Barke.2")


#Filter for genes present in c52_17 and c52_Barke and not in c17_Barke
interest <- c52_17 %>% filter(target %in% c52_Barke$target & !(target %in% c17_Barke$target))

interest_all <- inner_join(interest, BaRTv2, by='target')

write.csv(interest_all,paste0(result.folder,"/candidate_gene_complete_info50.csv"))
#Include information

#How many genes expressed in dataset?

length(unique(rownames(genes_3D_stat$DE.pval))) #19802


#DTU analysis - find DTU transcripts for 3H locus

DTUc52_17 <- DTU_trans %>% filter(contrast == "HEB_124_52.2-HEB_124_17.2")

DTUc52_Barke <- DTU_trans %>% filter(contrast == "HEB_124_52.2-Barke.2")

DTUc17_Barke <- DTU_trans %>% filter(contrast == "HEB_124_17.2-Barke.2")

#Filter for transcripts present in c52_17 and c52_Barke and not in c17_Barke
DTUinterest <- DTUc52_17 %>% filter(target %in% DTUc52_Barke$target & !(target %in% DTUc17_Barke$target))


BaRTv2_candidates <- read_csv(candidate_genes)


candidate_DTU <- DTUinterest %>% mutate(str_split_fixed(target,fixed("."), n = 2)[,1]) %>% rename(genes = 7) %>% filter(genes %in% BaRTv2_candidates$gene)#no DTU transcripts :(

#Try with DE transcripts

DEc52_17 <- DE_trans %>% filter(contrast == "HEB_124_52.2-HEB_124_17.2")

DEc52_Barke <- DE_trans %>% filter(contrast == "HEB_124_52.2-Barke.2")

DEc17_Barke <- DE_trans %>% filter(contrast == "HEB_124_17.2-Barke.2")

#Filter for transcripts present in c52_17 and c52_Barke and not in c17_Barke
DEinterest <- DEc52_17 %>% filter(target %in% DEc52_Barke$target & !(target %in% DEc17_Barke$target))


BaRTv2_candidates <- read_csv(candidate_genes)


candidate_DE <- DEinterest %>% mutate(str_split_fixed(target,fixed("."), n = 2)[,1]) %>% rename(genes = 6) %>% filter(genes %in% BaRTv2_candidates$gene)#

#DAS_genes

#Try and look at DAS genes:

DASc52_17 <- DAS_genes %>% filter(contrast == "HEB_124_52.2-HEB_124_17.2")

DASc52_Barke <- DAS_genes %>% filter(contrast == "HEB_124_52.2-Barke.2")

DASc17_Barke <- DAS_genes %>% filter(contrast == "HEB_124_17.2-Barke.2")

DASinterest <- DASc52_17 %>% filter(target %in% DASc52_Barke$target & !(target %in% DASc17_Barke$target))

BaRTv2_candidates <- read_csv(candidate_genes)


candidate_DAS <- DASinterest %>% filter(target %in% BaRTv2_candidates$gene)#

#None in region of interest

write.csv(DASinterest,paste0(result.folder,"/DAS_52-17_50subset.csv"))