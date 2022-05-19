#############################################################
#
# Ref to the ARTICLE 
# 
# Code to compute calculations presented in: https://www.biorxiv.org/content/10.1101/2021.12.20.472907v3
#
# Figure 4 and Supplementary Database 4
#  
# max.coulter@xelect.co.uk
# d.bulgarelli@dundee.ac.uk 
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#set working directory
setwd("/cluster/db/R_shared/3H_manuscript/")

#load required package function
library(tximport)
#This code uses the 3D RNA-seq app https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_command_line_user_manual.md
source("3D_source_code.R")
#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

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
input.folder <- 'salmon'
quant.folder <- 'quants'

quant_method <- 'salmon' # abundance generator
tximport_method <- 'lengthScaledTPM' # method to generate expression in tximport


################################################################################
#parameters for data pre-processing
################################################################################

#has sequencing technical replicates?
has_srep <- F

# parameter for low expression filters
cpm_cut <- 2
cpm_samples_n <- 3

# data normalisation parameter
norm_method <- 'TMM' ## norm_method is one of 'TMM','RLE' and 'upperquartile'


################################################################################
#parameters for DE gene analysis:
#EdgeR generalised linear model quasi-likelihood (glmQL)
#P-values were corrected using the Benjamini-Hochberg method to correct the false discovery rate. 
#Genes were considered to be DE if they had an adjusted P-value < 0.01 and a Log2FC >=1 or <=-1
################################################################################

pval_adj_method <- 'BH'
pval_cut <- 0.01
l2fc_cut <- 1
DE_pipeline = 'glmQL'

################################################################################
#Data genearation (I)
#metadata table with sample information
#used the same syntax for 3D app for simplicity
################################################################################

metatable <- read.csv('metadata-rhizosphere_RNA-Seqexperiment2only_DB.csv')
##select the columns of experimental design
factor_col <-'Tissue'
brep_col <- 'brep'
srep_col <- 'srep'
quant_col <- 'quant.files'
experiment <- 'experiment'
sample_date <- 'sampling.date'

#arrange information in the metatable
metatable$label <- as.vector(interaction(metatable[,factor_col]))
metatable$sample.name <- as.vector(interaction(metatable[,c(factor_col,brep_col,srep_col)]))
metatable$quant.folder <- file.path(quant.folder,metatable$quant.files, ifelse(quant_method=='salmon','quant.sf','abundance.h5'))


#Transcript-gene association mapping using BaRT2v18 as a reference
mapping <-read.csv("/cluster/db/mecoulter/RNAseq2/transcripts_genes_BaRT2v18.csv")
mapping <- data.frame(as.matrix(mapping),stringsAsFactors = F)
rownames(mapping) <- mapping$TXNAME

################################################################################
#Data genearation (III)
#merge sequencing replicates: no sequencing replicates, genes_counts and trans_counts remain the same by 
################################################################################

if(has_srep){
  idx <- paste0(metatable$label,'.',metatable[,brep_col])
  genes_counts <- sumarrays(genes_counts,group = idx)
  trans_counts <- sumarrays(trans_counts,group = idx)
  metatable_new <- metatable[metatable[,srep_col]==metatable[,srep_col][1],]
} else {
  metatable_new <- metatable
}

################################################################################
#Data genearation (III)
#generate transcript andgene expression
################################################################################

#Import transcript-level abundances and estimated counts for gene-level analysis packages

#transcript level
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

#gene level
txi_genes <- tximport(metatable$quant.folder,dropInfReps = T,
                      type = quant_method, tx2gene = mapping,
                      countsFromAbundance = tximport_method)

#assign names to columns of the dataset
colnames(txi_genes$counts) <-
  colnames(txi_genes$abundance) <-
  colnames(txi_genes$length) <-metatable$sample.name

## save the data
write.csv(txi_genes$counts,file=paste0(result.folder,'/counts_genes.csv'))
write.csv(txi_genes$abundance,file=paste0(result.folder,'/TPM_genes.csv'))
save(txi_genes,file=paste0(data.folder,'/txi_genes.RData'))

##extract gene and transcript read counts
genes_counts <- txi_genes$counts
trans_counts <- txi_trans$counts
trans_TPM <- txi_trans$abundance

################################################################################
#Data genearation (IV)
#Filter low-expression genes (see lines 65-73 for individual parameters)
################################################################################

#filtering at transcript level: A gene is expressed if any of its transcript is expressed
target_high <- low.expression.filter(abundance = trans_counts,
                                     mapping = mapping,
                                     abundance.cut = cpm_cut,
                                     sample.n = cpm_samples_n,
                                     unit = 'counts',
                                     Log=F)


#save expressed genes and transcripts
save(target_high,file=paste0(data.folder,'/target_high.RData'))

################################################################################
#Data genearation (V)
#Data normalisation: see line 73 for nomalisation method
################################################################################

#transcript level
dge <- DGEList(counts=trans_counts[target_high$trans_high,],
               group = metatable_new$label,
               genes = mapping[target_high$trans_high,])
trans_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(trans_dge,file=paste0(data.folder,'/trans_dge.RData'))

#gene level
dge <- DGEList(counts=genes_counts[target_high$genes_high,],
               group = metatable_new$label)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(genes_dge,file=paste0(data.folder,'/genes_dge.RData'))

################################################################################
#Data information
################################################################################

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

################################################################################
#Differential gene analysis
################################################################################

#set up the model
design <- model.matrix(~0 + Tissue + sampling.date, data=metatable_new)

#pair-wise contrast groups
contrast <- c('TissueHEB_124_17-TissueBarke','TissueHEB_124_52-TissueBarke','TissueHEB_124_52-TissueHEB_124_17')

#edgeR glmQL pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                design = design,
                                deltaPS = NULL,
                                contrast = contrast,
                                diffAS = F,
                                method = 'glmQL',
                                adjust.method = pval_adj_method)

#save the results
DDD.data$genes_3D_stat <- genes_3D_stat

#Summary DE genes
DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                            cutoff = c(adj.pval=pval_cut,
                                       log2FC=l2fc_cut))

DDD.data$DE_genes <- DE_genes

#save the results
write.csv(DE_genes,file=paste0(result.folder,'/DE genes.csv'),row.names = F)

################################################################################
#Figure 4 generation
################################################################################

#extract genes differentially regulated in 124_17-Barke comparison
Bvs17 <- DE_genes[(rownames(DE_genes)[which(DE_genes$contrast == 'TissueHEB_124_17-TissueBarke')]), ]
#identify the genes
Bvs17_genes <- Bvs17$target 
#number of genes
length(Bvs17_genes)
#no gene passing the threshold criteria

#extract genes differentially regulated in 124_52-Barke comparison
Bvs52 <- DE_genes[(rownames(DE_genes)[which(DE_genes$contrast == 'TissueHEB_124_52-TissueBarke')]), ]
#identify the genes
Bvs52_genes <- Bvs52$target 
#number of genes
length(Bvs52_genes)

#extract genes differentially regulated in 124_52-124_17 comparison
L17vs52 <- DE_genes[(rownames(DE_genes)[which(DE_genes$contrast == 'TissueHEB_124_52-TissueHEB_124_17')]), ]
#identify the genes
L17vs52_genes <- L17vs52$target 
#number of genes
length(L17vs52_genes)

#Unique genes Barke vs 124_52
Bvs52_genes_unique <- setdiff(Bvs52_genes, L17vs52_genes)
length(Bvs52_genes_unique)
#50

#Unique genes 124_17 vs 124_52
L17vs52_genes_unique <- setdiff(L17vs52_genes, Bvs52_genes)
length(L17vs52_genes_unique)
#3

#intersection
L52_vs_elite <- intersect(L17vs52_genes, Bvs52_genes)
length(L52_vs_elite)
#34

#venn diagram created in illustrator

################################################################################
#Supplementary database 4 generation
################################################################################

#subset the 34 genes from L17vs52
candidate_info <- L17vs52[, 2:5]
rownames(candidate_info) <- L17vs52$target
candidate_info <- candidate_info[L52_vs_elite, ]

#upload gene annotation
gene_annotation <- read.delim("/cluster/db/mecoulter/BaRT2v18/BaRT_2_18_annotation_genes.txt", header = T, row.names = 1)
colnames(gene_annotation)

#create the dataset
gene_annotation_candidates <-gene_annotation[rownames(candidate_info), ] 
gene_annotation_candidates[1:5, ]
gene_annotation_candidates_2 <- gene_annotation_candidates[, 1:4]
colnames(gene_annotation_candidates_2) <- c("Chromosome", "Start", "End", "Strand")
gene_annotation_candidates_panzer <- as.data.frame(gene_annotation_candidates[, "Pannzer.annotation"])
colnames(gene_annotation_candidates_panzer) <- c("Pannzer annotation")
rownames(gene_annotation_candidates_panzer) <- row.names(gene_annotation_candidates_2) 
gene_annotation_candidates_panzer

gene_annotation_candidates_info <- cbind(candidate_info, gene_annotation_candidates_2)
gene_annotation_candidates_info_pannzer <- cbind(gene_annotation_candidates_info, gene_annotation_candidates_panzer)
gene_annotation_candidates_info_pannzer[1:5, ]
#save the file to generate Supplementary Database 4

#end
