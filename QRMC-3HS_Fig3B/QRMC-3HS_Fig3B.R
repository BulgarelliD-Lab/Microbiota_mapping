#####################################################################################
#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:https://www.biorxiv.org/content/10.1101/2021.12.20.472907v1
#  Figure 3B 
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

library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("dplyr")
library("scales")
library("VennDiagram")
library("tidyverse")
library("forcats")
library("PMCMRplus")
library("phyloseq")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

#################################################
################################################
##Fig. 3B  Fungal community
###############################################
##############################################

#Import data picked with  JH23 DADA2
JH23_data_phyloseq_DADA2 <-readRDS("JH23_dada2.rds")
JH23_data_phyloseq_sample <- JH23_data_phyloseq_DADA2
JH23_data_phyloseq_sample

sample_names(JH23_data_phyloseq_sample)
rank_names(JH23_data_phyloseq_sample)
JH23_tax_table<-tax_table(JH23_data_phyloseq_sample)

dim(tax_table(JH23_data_phyloseq_sample))

##################################################################
#Pre-processing: remove Chloroplast and Mitochondria but retain NA
#################################################################
JH23_taxa_raw<-tax_table(JH23_data_phyloseq_sample)
#write.table(JH23_taxa_raw, file="JH23_taxa_raw_taxonomy_info.txt", sep="\t")

#Remove chloroplasts but retain"NA" JH23
JH23_no_chlor <-subset_taxa(JH23_data_phyloseq_DADA2, (Order!="Chloroplast") | is.na(Order))
JH23_no_chlor

#Remove mitochondria but retains "NA"
JH23_no_plants <-subset_taxa(JH23_no_chlor, (Family!="Mitochondria") | is.na(Family))
JH23_no_plants


colnames (tax_table(JH23_no_plants))
rownames (tax_table(JH23_no_plants))
dim(tax_table(JH23_no_plants))

#########################################################################
# Remove ASVs assigned to NA at phylum level
#########################################################################

JH23_phyloseq_DADA2_no_plants_1 <- subset_taxa (JH23_no_plants, Phylum!= "NA")
JH23_phyloseq_DADA2_no_plants_1

#########################################################################
# Abundance threshold: 20 reads, set as 2% the minimum number of samples
#########################################################################

JH23_phyloseq_DADA2_no_plants_10K <- prune_samples(sample_sums(JH23_phyloseq_DADA2_no_plants_1)>=10000, JH23_phyloseq_DADA2_no_plants_1)
JH23_phyloseq_DADA2_no_plants_10K


JH23_phyloseq_DADA2_no_plants_2 = filter_taxa(JH23_phyloseq_DADA2_no_plants_10K, function(x) sum(x > 20) > (0.02 *length(x)), TRUE)
JH23_phyloseq_DADA2_no_plants_2 
sort(sample_sums(JH23_phyloseq_DADA2_no_plants_2))

JH23_taxa<-tax_table(JH23_phyloseq_DADA2_no_plants_2 )
#write.table(JH23_taxa, file="JH23_taxa_info.txt", sep="\t")

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH23_phyloseq_DADA2_no_plants_2))/sum(sample_sums(JH23_phyloseq_DADA2_no_plants_1))*100
ratio

#########################################################################
# Import a simplified mapping file
########################################################################

#Import the simplified mapping file
design <- read.delim("JH23_Map.txt", sep = "\t", header=TRUE, row.names=1)

JH23_map <- sample_data(design) 

sample_data(JH23_phyloseq_DADA2_no_plants_2)<-JH23_map
JH23_phyloseq_DADA2_no_plants_2
#########################################################################
# Aggregate samples at genus and family level (Note: NArm set to false to keep NA at Genus level)
########################################################################

JH23_data_phyloseq_genus <- tax_glom(JH23_phyloseq_DADA2_no_plants_2, taxrank= "Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

JH23_genus_taxa<-tax_table(JH23_data_phyloseq_genus)
#write.table(JH23_genus_taxa, file="JH23_genus_taxa_taxonomy_info.txt", sep="\t")

JH23_data_phyloseq_family <- tax_glom(JH23_phyloseq_DADA2_no_plants_2, taxrank= "Family", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

#Compare the two objects
#ASVs
JH23_phyloseq_DADA2_no_plants_2
sort(sample_sums(JH23_phyloseq_DADA2_no_plants_2))


#Genus
JH23_data_phyloseq_genus
sort(sample_sums(JH23_data_phyloseq_genus))


#Family
JH23_data_phyloseq_family
sort(sample_sums(JH23_data_phyloseq_family))

####################################################################################
# Rarefy at an even sequencing depth (10,600) and "freeze" these objects for downstream analyses
#####################################################################################

#ASVs : ignore the warnings, the object will be saved right after
JH23_rare_ASV_10K <- rarefy_even_depth(JH23_phyloseq_DADA2_no_plants_2, 10200)
#saveRDS(JH23_rare_ASV_10K, file ="JH23_rare_ASV_10K.rds")

JH23_rare_ASV_10K<- readRDS("JH23_rare_ASV_10K.rds")


######################################################################
######################################################################
## Figure 3B_Constrained ordination
#####################################################################
#####################################################################

JH23_rhizosphere<- subset_samples(JH23_rare_ASV_10K, Genotype == "Barke" |Genotype == "52"|Genotype == "17") 
JH23_rhizosphere

JH23_rare_sqrt<-transform_sample_counts(JH23_rhizosphere, function (x) sqrt(x)) 


#constrained ordination: constrained for genotype
#rarefied data
JH23_rhizosphere_CAP <- ordinate(JH23_rare_sqrt, "CAP", "bray", ~ Genotype)
plot_ordination(JH23_rare_sqrt, JH23_rhizosphere_CAP, color = "Genotype")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH23_rare_sqrt, JH23_rhizosphere_CAP , color = "Genotype")
p = p + geom_point(size = 6, alpha = 0.80, shape =17)
p = p + scale_colour_manual(values = c("#56B4E9","#E69F00","#0072B2"))
p + ggtitle("CAP 16S data, Bray distance-rarefied samples")

#extract the mapping file
design_rhizosphere <- design[colnames(otu_table(JH23_rhizosphere)), ]

#BC distance
BC <- phyloseq::distance(JH23_rhizosphere, "bray")

adonis(BC ~ Genotype, data= design_rhizosphere , permutations = 5000)

adonis(BC ~ Batch*Genotype , data=design_rhizosphere , permutation =5000)

############################################################################################################################
##End