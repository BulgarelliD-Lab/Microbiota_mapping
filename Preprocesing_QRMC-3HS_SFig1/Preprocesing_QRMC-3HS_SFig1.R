#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:https://www.biorxiv.org/content/10.1101/2021.12.20.472907v3
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
library("vegan")
library ("ape")
library("PMCMR")
library("plyr")
library("grid")


#retrieve R and package versions for the reproducibility of the code
sessionInfo()

#set the working directory
setwd()
getwd()

#import data picked with JH07 DADA2
JH07_data_phyloseq_DADA2 <-readRDS("JH07_dada2.silva138.rds")
JH07_map <- read.delim("Map_JH07.txt", row.names = 1)
sample_data(JH07_data_phyloseq_DADA2) <- JH07_map
JH07_data_phyloseq_sample <- subset_samples(JH07_data_phyloseq_DADA2, Type =="Sample")
JH07_data_phyloseq_sample


rank_names(JH07_data_phyloseq_sample)
tax_table(JH07_data_phyloseq_sample)[1:1, 1:6]
dim(tax_table(JH07_data_phyloseq_sample))


#import data picked with JH15 DADA2
JH15_data_phyloseq_DADA2 <-readRDS("JH15_dada2.silva138.rds")
JH15_data_phyloseq_sample <- subset_samples(JH15_data_phyloseq_DADA2)
JH15_data_phyloseq_sample


sample_names(JH15_data_phyloseq_sample)
rank_names(JH15_data_phyloseq_sample)
tax_table(JH15_data_phyloseq_sample)[1:1, 1:6]
dim(tax_table(JH15_data_phyloseq_sample))

#Pre-processing: remove Chloroplast and Mitochondria but retain NA
####################################################################
#remove chloroplasts but retain"NA" JH07
JH07_no_chlor <-subset_taxa(JH07_data_phyloseq_sample, (Order!="Chloroplast") | is.na(Order))
JH07_no_chlor

#remove mitochondria but retains "NA"
JH07_no_plants <-subset_taxa(JH07_no_chlor, (Family!="Mitochondria") | is.na(Family))
JH07_no_plants

colnames (tax_table(JH07_no_plants))
rownames (tax_table(JH07_no_plants))
dim(tax_table(JH07_no_plants))

#remove chloroplasts but retain"NA" JH15
JH15_no_chlor <-subset_taxa(JH15_data_phyloseq_DADA2, (Order!="Chloroplast") | is.na(Order))
JH15_no_chlor

#remove mitochondria but retains "NA"
JH15_no_plants <-subset_taxa(JH15_no_chlor, (Family!="Mitochondria") | is.na(Family))
JH15_no_plants


colnames (tax_table(JH15_no_plants))
rownames (tax_table(JH15_no_plants))
dim(tax_table(JH15_no_plants))

#Merging Phyloseq objects
########################################################
JH15_JH07_phloseq_DADA2_no_plants <- merge_phyloseq(JH15_no_plants, JH07_no_plants)
JH15_JH07_phloseq_DADA2_no_plants

#Prune putative contaminant ASVs
#######################################################
#Import the list of contaminant ASVs from JH06 library
Contaminant_ASVs <- read.delim("JH06_contaminant_list_2.txt", header = FALSE)

#Identify the proportion of putative contaminants in the merged object
JH15_JH07_contaminants <- intersect(taxa_names(JH15_JH07_phloseq_DADA2_no_plants),Contaminant_ASVs)
JH15_JH07_contaminants

#Remove ASVs assigned to NA at phylum level
#######################################################
JH15_JH07_phloseq_DADA2_no_plants_1 <- subset_taxa(JH15_JH07_phloseq_DADA2_no_plants, Phylum!= "NA")
JH15_JH07_phloseq_DADA2_no_plants_1


#Abundance threshold: 20 reads, set as 2% the minimum number of samples
########################################################################
JH15_JH07_phloseq_DADA2_no_plants_10K <- prune_samples(sample_sums(JH15_JH07_phloseq_DADA2_no_plants_1)>=10000, JH15_JH07_phloseq_DADA2_no_plants_1)
JH15_JH07_phloseq_DADA2_no_plants_10K


JH15_JH07_phloseq_DADA2_no_plants_2 = filter_taxa(JH15_JH07_phloseq_DADA2_no_plants_10K, function(x) sum(x > 20) > (0.02 *length(x)), TRUE)
JH15_JH07_phloseq_DADA2_no_plants_2 
sort(sample_sums(JH15_JH07_phloseq_DADA2_no_plants_2))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH15_JH07_phloseq_DADA2_no_plants_2))/sum(sample_sums(JH15_JH07_phloseq_DADA2_no_plants_1))*100
ratio

#Generate a simplified mapping file
#####################################

#Import the simplified mapping file
JH15_JH07_map <- read.delim("JH15_JH07_Map.txt", row.names = 1)


#and replace the old mapping file
sample_data(JH15_JH07_phloseq_DADA2_no_plants_2) <- JH15_JH07_map 
sample_data(JH15_JH07_phloseq_DADA2_no_plants_2)


#Aggregate samples at genus and family level (Note: NArm set to false to keep NA at Genus level)
#############################################
JH15_JH07_data_phyloseq_genus <- tax_glom(JH15_JH07_phloseq_DADA2_no_plants_2, taxrank= "Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
JH15_JH07_data_phyloseq_family <- tax_glom(JH15_JH07_phloseq_DADA2_no_plants_2, taxrank= "Family", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

#Compare the objects
####################
#ASVs
JH15_JH07_phloseq_DADA2_no_plants_2
sort(sample_sums(JH15_JH07_phloseq_DADA2_no_plants_2))
#Genus
JH15_JH07_data_phyloseq_genus
sort(sample_sums(JH15_JH07_data_phyloseq_genus))
#Family
JH15_JH07_data_phyloseq_family
sort(sample_sums(JH15_JH07_data_phyloseq_family))

# Rarefy at an even sequencing depth (10,600) and "freeze" these objects for downstream analyses
#####################################################################################

#ASVs : ignore the warnings, the object will be saved right after
JH15_JH07_rare_10K <- rarefy_even_depth(JH15_JH07_phloseq_DADA2_no_plants_2, 10600)
saveRDS(JH15_JH07_rare_10K, file ="JH15_JH07_rare_ASV_10K.rds")

#JH15_JH07_rare_10K<- readRDS("JH15_JH07_rare_10K.rds")

#Genus : ignore the warnings, the object will be saved right after
JH15_JH07_rare_genus_10K <- rarefy_even_depth(JH15_JH07_data_phyloseq_genus, 10600)
saveRDS(JH15_JH07_rare_genus_10K, file ="JH15_JH07_rare_genus_10K.rds")

#JH15_JH07_rare_genus_10K<- readRDS("JH15_JH07_rare_genus_10K.rds")

#Family : ignore the warnings, the object will be saved right after
JH15_JH07_rare_family_10K <- rarefy_even_depth(JH15_JH07_data_phyloseq_family, 10600)
saveRDS(JH15_JH07_rare_family_10K, file ="JH15_JH07_rare_family_10K.rds")

#JH15_JH07_rare_family_10K<- readRDS("JH15_JH07_rare_family_10K.rds")

####End

######################################################################################

##Supplementary_Figure 1 Tecnical reproducibility

######################################################################################

#Subset for the control samples for the determination of the core measurable microbiome
JH07_data_phyloseq_CMM <- subset_samples(JH07_data_phyloseq_DADA2, CMM =="Yes")
JH07_data_phyloseq_CMM
JH07_data_phyloseq_sample <- subset_samples(JH07_no_plants, Type =="Sample")
JH07_data_phyloseq_sample

#How to integrate the tecnical reproducibility into the total counts. 
#Using rarefied counts for technical reproducibility is problematic, unless we use different thresholds!
#e.g., 10K, 20K, 50K, the lowest common point will be used to set the reproducibility of the entire dataset!

#Define the number of reads of the technical replicates
sample_sums(JH07_data_phyloseq_CMM)

#Rarefy at 10K, 25K and 50K
#10K
JH07_data_phyloseq_CMM_10K <- rarefy_even_depth(JH07_data_phyloseq_CMM, sample.size = 10000)
#25K
JH07_data_phyloseq_CMM_25K <- rarefy_even_depth(JH07_data_phyloseq_CMM, sample.size = 25000)
#50K
JH07_data_phyloseq_CMM_50K <- rarefy_even_depth(JH07_data_phyloseq_CMM, sample.size = 50000)

#Extract the data
JH07_data_phyloseq_CMM_10K_table <- as.data.frame(otu_table(JH07_data_phyloseq_CMM_10K))
JH07_data_phyloseq_CMM_25K_table <- as.data.frame(otu_table(JH07_data_phyloseq_CMM_25K))
JH07_data_phyloseq_CMM_50K_table <- as.data.frame(otu_table(JH07_data_phyloseq_CMM_50K))

#save the file for the reproducibility of the code
#write.table(JH07_data_phyloseq_CMM_10K_table, file="JH07_data_phyloseq_CMM_10K_table.txt", sep="\t")
#write.table(JH07_data_phyloseq_CMM_25K_table, file="JH07_data_phyloseq_CMM_25K_table.txt", sep="\t")
#write.table(JH07_data_phyloseq_CMM_50K_table, file="JH07_data_phyloseq_CMM_50K_table.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_CMM_10K <- read.delim("JH07_data_phyloseq_CMM_10K_table2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
dat_count_CMM_25K <- read.delim("JH07_data_phyloseq_CMM_25K_table2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
dat_count_CMM_50K <- read.delim("JH07_data_phyloseq_CMM_50K_table2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)


#calculate correlations
pw_cor <- cbind (cval12, cval13,cval32)
# 5 reads
apply(pw_cor,1,mean)[5]
# 10 reads#make different thresholds to calculate library sizes
ts <- (1:50)
cval12<-c()
#Calculate correlation at different increasing abundance for each pair of technical replicate
#10K
for(threshold in ts){
  #select OTUs with relative abundance higher than threshold
  datr <- dat_count_CMM_10K[apply(dat_count_CMM_10K[,c(1,2)],1,max) > threshold, ]
  cval12<- c(cval12,cor(datr[,1],datr[,2],meth = "spearman"))
}
cval13<-c()
for(threshold in ts){
  datr <- dat_count_CMM_10K[apply(dat_count_CMM_10K[,c(1,3)],1,max) > threshold, ]
  cval13<- c(cval13,cor(datr[,1],datr[,3],meth = "spearman"))
}
cval32<-c()
for(threshold in ts){
  datr <- dat_count_CMM_10K[apply(dat_count_CMM_10K[,c(3,2)],1,max) > threshold, ]
  cval32<- c(cval32,cor(datr[,3],datr[,2],meth = "spearman"))
}

#set the screen
dev.off()
#par( mfrow = c( 3, 1 ) )
plot(cval12,x=ts, xlim =c(0,50), yaxp=c(-1,1,20),xlab="ASV minimal abundance (reads)",ylab=c("Spearman rank correlation"), main = "ASV technical reproducibility sample size 10K")
points(cval13,x=ts,col="#D55E00",pch=19)
points(cval32,x=ts,col="#009E73",pch=19)
points(cval12,x=ts,col="black", pch=19)
#calculate correlations
pw_cor <- cbind (cval12, cval13,cval32)
# 5 reads
apply(pw_cor,1,mean)[5]
# 10 reads
apply(pw_cor,1,mean)[10]
#15 reads
apply(pw_cor,1,mean)[15]
#20 reads
apply(pw_cor,1,mean)[20]

#Calculate correlation at different increasing abundance for each pair of technical replicate
#25K
cval12<-c()
for(threshold in ts){
  #select OTUs with relative abundance higher than threshold
  datr <- dat_count_CMM_25K[apply(dat_count_CMM_25K[,c(1,2)],1,max) > threshold, ]
  cval12<- c(cval12,cor(datr[,1],datr[,2],meth = "spearman"))
}
cval13<-c()
for(threshold in ts){
  datr <- dat_count_CMM_25K[apply(dat_count_CMM_25K[,c(1,3)],1,max) > threshold, ]
  cval13<- c(cval13,cor(datr[,1],datr[,3],meth = "spearman"))
}
cval32<-c()
for(threshold in ts){
  datr <- dat_count_CMM_25K[apply(dat_count_CMM_25K[,c(3,2)],1,max) > threshold, ]
  cval32<- c(cval32,cor(datr[,3],datr[,2],meth = "spearman"))
}
plot(cval12,x=ts, xlim =c(0,50), yaxp=c(-1,1,20),xlab=" ASV minimal abundance (reads)",ylab=c("Spearman rank correlation"), main = "ASV technical reproducibility sample size 25K" )
points(cval13,x=ts,col="red")
points(cval32,x=ts,col="blue")


#calculate correlations
pw_cor <- cbind (cval12, cval13,cval32)
#5 reads
apply(pw_cor,1,mean)[10]
#10 reads
apply(pw_cor,1,mean)[10]
#15 reads
apply(pw_cor,1,mean)[15]
#20 reads
apply(pw_cor,1,mean)[20]
#30 reads
apply(pw_cor,1,mean)[30]

#Calculate correlation at different increasing abundance for each pair of technical replicate
#50K
cval12<-c()
for(threshold in ts){
  #select OTUs with relative abundance higher than threshold
  datr <- dat_count_CMM_50K[apply(dat_count_CMM_50K[,c(1,2)],1,max) > threshold, ]
  cval12<- c(cval12,cor(datr[,1],datr[,2],meth = "spearman"))
}
cval13<-c()
for(threshold in ts){
  datr <- dat_count_CMM_50K[apply(dat_count_CMM_50K[,c(1,3)],1,max) > threshold, ]
  cval13<- c(cval13,cor(datr[,1],datr[,3],meth = "spearman"))
}
cval32<-c()
for(threshold in ts){
  datr <- dat_count_CMM_50K[apply(dat_count_CMM_50K[,c(3,2)],1,max) > threshold, ]
  cval32<- c(cval32,cor(datr[,3],datr[,2],meth = "spearman"))
}
plot(cval12,x=ts, xlim =c(0,50), yaxp=c(-1,1,20),xlab=" ASV minimal abundance (reads)",ylab=c("Spearman rank correlation"), main = "ASV technical reproducibility sample size 50K" )
points(cval13,x=ts,col="red")
points(cval32,x=ts,col="blue")

#calculate correlations
pw_cor <- cbind (cval12, cval13,cval32)
# 5 reads
apply(pw_cor,1,mean)[5]
# 10 reads
apply(pw_cor,1,mean)[10]
#15 reads
apply(pw_cor,1,mean)[15]
#20  reads
apply(pw_cor,1,mean)[25]
#30 reads
apply(pw_cor,1,mean)[30]

####End
###########################################################################################################
