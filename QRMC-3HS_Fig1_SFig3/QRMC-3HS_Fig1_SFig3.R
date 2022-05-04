#####################################################################################
#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:https://www.biorxiv.org/content/10.1101/2021.12.20.472907v1
#  Figure 1, S Figure 3 
#   
#  c.m.z.escuderomartinez@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#required packages 
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("dplyr")
library("grid")
library("magrittr")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

##################################################################################
## Figure 1 
##################################################################################
################################################################################

#Import the file: 
JH07_JH15_ASV <- readRDS("JH15_JH07_rare_ASV_10K.rds")
JH07_JH15_ASV


#Subset the library of interest, in this case JH07
JH07_ASV <- subset_samples(JH07_JH15_ASV, Library =="JH07")
JH07_ASV

JH07_no_plants_dada2<-tax_table(JH07_ASV)


#Assign the new mapping
########################
#import the simplified mapping file
design <- read.delim("Map_JH07_phyloseq_3.txt", row.names = 1)


#and replace the old mapping file
sample_data(JH07_ASV) <- design 
sample_data(JH07_ASV)

################################################
################################################
## Figure 1A_Ternary plot
################################################
################################################


#Subset the data for the parental lines and the parental rhizospheres
######################################################################
JH07_parentals <- subset_samples(JH07_ASV, Description != "BC1S3")

JH07_parentals_rhizo <- subset_samples(JH07_ASV, Description == "Barke"| Description == "HID144")


#Calculation of differential abundant taxa between the parental lines using Deseq2
######################################################################################

#Extract count data 
JH07_counts <- otu_table(JH07_parentals)
countData = as.data.frame(JH07_counts)
colnames(JH07_counts)

#The design file containing sample information
colData = design[colnames(JH07_counts), ]
rownames(colData)

#Construct a DESeq data set combining count data and sample information
JH07_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#Execute the differential count analysis with the function DESeq 
JH07_cds_test <- DESeq(JH07_cds, fitType="local", betaPrior = FALSE)

#saveRDS(JH07_cds_test, file="JH07_cds_test_ASVs")
JH07_cds_test<-readRDS(file="JH07_cds_test_ASVs")
levels(JH07_cds_test$Description)


#Define the OTUs significantly enriched in the rhizosphere samples
rhizo_Barke <- results(JH07_cds_test, contrast = c("Description", "Barke", "Bulk")) 
rhizo_HID144 <- results(JH07_cds_test, contrast = c("Description", "HID144", "Bulk")) 

#inspect a result file
rhizo_Barke  
mcols(rhizo_Barke  , use.names=TRUE)
rhizo_HID114  
mcols(rhizo_HID144  , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_Barke_FDR_005 <- rhizo_Barke[(rownames(rhizo_Barke)[which(rhizo_Barke$padj <0.05)]), ]
rhizo_HID144_FDR_005 <- rhizo_HID144[(rownames(rhizo_HID144)[which(rhizo_HID144$padj <0.05)]), ]

#enriched in each of the parental lines (positive fold change)
rhizo_Barke_enriched <-  rhizo_Barke[(rownames(rhizo_Barke)[which(rhizo_Barke$log2FoldChange > 0)]), ]
rhizo_HID144_enriched <-  rhizo_HID144[(rownames(rhizo_HID144)[which(rhizo_HID144$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_Barke_enriched_FDR005 <- intersect(rownames(rhizo_Barke_FDR_005), rownames(rhizo_Barke_enriched))
rhizo_HID144_enriched_FDR005 <- intersect(rownames(rhizo_HID144_FDR_005), rownames(rhizo_HID144_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_Barke_enriched_FDR005)
length(rhizo_HID144_enriched_FDR005)

#pair-wise comparison
Barke_HID144 <- results(JH07_cds_test, contrast = c("Description", "Barke", "HID144")) 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Barke_HID144_FDR_005 <- Barke_HID144[(rownames(Barke_HID144)[which(Barke_HID144$padj <0.05)]), ]

#enriched in barke
Barke_enriched_HID144 <-  Barke_HID144_FDR_005[(rownames(Barke_HID144_FDR_005)[which(Barke_HID144_FDR_005$log2FoldChange > 0)]), ]
nrow(Barke_enriched_HID144)

#enriched in HID144
HID144_enriched_Barke <-  Barke_HID144_FDR_005[(rownames(Barke_HID144_FDR_005)[which(Barke_HID144_FDR_005$log2FoldChange < 0)]), ]
nrow(HID144_enriched_Barke)

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Barke_enriched_HID144_FDR005 <- intersect(rownames(Barke_HID144_FDR_005), rownames(Barke_enriched_HID144))
HID144_enriched_Barke_FDR005 <- intersect(rownames(Barke_HID144_FDR_005), rownames(HID144_enriched_Barke))

#are these ASV Barke enriched from soil?
Barke_enriched_FDR005_vs_HID144_rhizo <- intersect(rhizo_Barke_enriched_FDR005, Barke_enriched_HID144_FDR005)
length(Barke_enriched_FDR005_vs_HID144_rhizo)

JH07_no_plants_dada2[Barke_enriched_FDR005_vs_HID144_rhizo, ]

#are these ASV HID144 enriched from soil?
HID144_enriched_FDR005_vs_Barke_rhizo <- intersect(rhizo_HID144_enriched_FDR005, HID144_enriched_Barke_FDR005)
length(HID144_enriched_FDR005_vs_Barke_rhizo)

JH07_no_plants_dada2[HID144_enriched_FDR005_vs_Barke_rhizo, ]

#retrieve taxonomy information
JH07_HID144_enriched<-JH07_no_plants_dada2[HID144_enriched_FDR005_vs_Barke_rhizo, ]
JH07_Barke_enriched<-JH07_no_plants_dada2[Barke_enriched_FDR005_vs_HID144_rhizo, ]

#write.table(JH07_HID144_ASVs_enriched, file="JH07_merged_HID144_ASV_enriched_taxa_dada2_v1220.txt", sep="\t")
#write.table(JH07_Barke_ASVs_enriched, file="JH07_merged_Barke_ASV_enriched_taxa_dada2_v1220.txt", sep="\t")

#####Calculation proportion of reads enriched in each of the parental lines respect the total rhizosphere reads
###median values
##########################################################################################################
#Remove bulk soil whose abundances are markedly different from plant-associated ones
JH07_plants <- subset_samples(JH07_ASV, Description !="Bulk")

###Proportion of reads in HID144

JH07_plants_HID144 <- prune_taxa(HID144_enriched_FDR005_vs_Barke_rhizo, JH07_plants)
JH07_plants_HID144

HID144_proportion <- median(sample_sums(JH07_plants_HID144))/median(sample_sums(JH07_plants)) * 100
HID144_proportion
###Proportion of reads in Barke

JH07_plants_Barke <- prune_taxa(Barke_enriched_FDR005_vs_HID144_rhizo, JH07_plants)
JH07_plants_Barke

Barke_proportion <- median(sample_sums(JH07_plants_Barke))/median(sample_sums(JH07_plants)) * 100
Barke_proportion

#Visualisation with ternary plot
#################################
#extract the means of Modern and Wild
JH07_base_mean <- sapply(levels(JH07_cds_test$Description), function(lvl) rowMeans(counts(JH07_cds_test,normalized=TRUE)[,JH07_cds_test$Description == lvl] ) )
head(JH07_base_mean)
### Ternary plot
mean_soil <- as.data.frame(JH07_base_mean[, 2])
mean_barke <- as.data.frame(JH07_base_mean[, 1])
mean_HID144 <- as.data.frame(JH07_base_mean[, 3])
temat_1A <- cbind(mean_HID144, mean_barke, mean_soil)
colnames(temat_1A) <- c("","","")#here the vertice names
dim(temat_1A)
# remove OTUs which are only detected in 
temat_1A <- temat_1A[!(rowSums(temat_1A) == 0),]
dim(temat_1A)
#colors and plotting
dev.off()
fig_colors <- ifelse(rownames(temat_1A) %in% Barke_enriched_FDR005_vs_HID144_rhizo, "#0072B2","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[HID144_enriched_FDR005_vs_Barke_rhizo] <- "#D55E00"
tern_e(temat_1A, prop=T, col=fig_colors, grid_color="black", labels_color="black", pch=19, main="")

#Count information
temat_1A[HID144_enriched_FDR005_vs_Barke_rhizo,]
temat_1A[Barke_enriched_FDR005_vs_HID144_rhizo, ]

#info
length(HID144_enriched_FDR005_vs_Barke_rhizo)
length(Barke_enriched_FDR005_vs_HID144_rhizo)


############################################
###########################################
## Figure 1B_Constrained ordination
###########################################
###########################################

#Constrained ordination: constrained for Description
#Rarefied data
JH07_CAP <- ordinate(JH07_ASV, "CAP", "bray", ~ Description)
plot_ordination(JH07_ASV, JH07_CAP, color = "Description")

#Assign shapes to Soil type and color to Sample type
p=plot_ordination(JH07_ASV, JH07_CAP , color = "Description")
p = p + geom_point(size = 5, alpha = 0.80)
p = p + scale_colour_manual(values = c("#0072B2", "#CC79A7", "#000000","#D55E00" ))
p + ggtitle("CAP 16S data, Bray distance-rarefied samples")

#ANOVA on the axis
anova(JH07_CAP, permutations=5000)

#Genotype effect (rhizosphere only)
JH07_rhizo <- subset_samples(JH07_ASV, Microhabitat == "Rhizosphere")
design_rhizosphere <- design[colnames(otu_table(JH07_rhizo)), ]
#BC distance
BC <- phyloseq::distance(JH07_rhizo, "bray")
adonis(BC ~ Description, data= design_rhizosphere, permutations = 5000)


################################################################################################
#Supplementary figure 3 Top 5 Phyla of differently enriched taxa between the parental lines
###############################################################################################

JH07_plants_HID144 <- prune_taxa(HID144_enriched_FDR005_vs_Barke_rhizo, JH07_ASV)
JH07_plants_HID144

JH07_plants_Barke <- prune_taxa(Barke_enriched_FDR005_vs_HID144_rhizo, JH07_ASV)
JH07_plants_Barke

#Merge the two Phyloseq objects with differential enriched (DE) between the parental lines
DE_parentals_parentals_rhizo<-merge_phyloseq(JH07_plants_HID144,JH07_plants_Barke)
DE_parentals_parentals_rhizo_2<-subset_samples(DE_parentals_parentals_rhizo, Description == "Barke"| Description == "HID144")
## agglomerate at Phylum taxonomic rank, this Is a new phyloseq object
DE_parentals_phylum <- (x1 <- tax_glom(DE_parentals_parentals_rhizo_2, taxrank="Phylum") )
# ## How many taxa before/after agglomeration?
ntaxa(DE_parentals_phylum); ntaxa(x1)

tax_table(DE_parentals_phylum)[1:7, ]

sample_data(DE_parentals_phylum)
#transform in to relative abundance
ps_phylum_0 <- transform_sample_counts(DE_parentals_phylum, function(x) x / sum(x))
#abundance of all samples plot
plot_bar(ps_phylum_0, fill="Phylum")
#merge samples by treatment
ps_phylum_1 <- merge_samples(ps_phylum_0, "Genotype")
#transform to relative abudance
ps_phylum_2 <- transform_sample_counts(ps_phylum_1, function(x) x / sum(x))
#plot ra of all phyla
plot_bar(ps_phylum_2, fill="Phylum")
#plotting based on sample type
df_phylum <- psmelt(ps_phylum_2)
#write.table(df_phylum, file="df_phylum.txt ", sep="\t")
sample_data(ps_phylum_0)
tax_table(ps_phylum_2)


plot_bar(ps_phylum_2, fill="Phylum")

sample_data(ps_phylum_2)
tax_table(ps_phylum_2)

#plotting based on sample type
df_phylum <- psmelt(ps_phylum_2)
#write.table(df_phylum, file="df_phylum.txt ", sep="\t")


top_phylum <- df_phylum
  group_by("Description", Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phylum

#write.table(top_phylum, file="top_phylum.txt ", sep="\t")
top5_phylum <- top_phylum$Phylum[1:5]
df_phylum_0 <- df_phylum %>%
  mutate(Phylum = fct_other(Phylum, top5_phylum))
plot_top_5_phylum <-ggplot(df_phylum_0, aes(Sample, Abundance, fill = fct_reorder(Phylum, Abundance))) + geom_col() + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
plot_top_5_phylum + scale_fill_manual(values = c("#999999","#CC79A7","#0072B2","#F0E442","#009E73","#D55E00"))


############################################################################################################################
##End