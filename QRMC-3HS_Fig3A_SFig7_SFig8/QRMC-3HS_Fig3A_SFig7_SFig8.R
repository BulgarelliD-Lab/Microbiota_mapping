#####################################################################################
#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in: https://www.biorxiv.org/content/10.1101/2021.12.20.472907v3
#  Figure 3A, SFigure 7AD,  SFigure 8 qPCR
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
library("tidyverse")
library("forcats")
library("PMCMRplus")
library("phyloseq")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

################################################################################
##################################################################################
## Figure 3A
##################################################################################
##################################################################################
#import the file: 
JH07_JH15_ASV <- readRDS("JH15_JH07_rare_ASV_10K.rds")
JH07_JH15_ASV


#Subset the library of interest, in this case JH15
JH15_ASV <- subset_samples(JH07_JH15_ASV, Library =="JH15")
JH15_ASV

JH15_no_plants_dada2<-tax_table(JH15_ASV)

#########################################################################
#Assign the new mapping
########################################################################

#import the simplified mapping file
design <- read.delim("JH15_Map_batch.txt", row.names = 1)


#and replace the old mapping file
sample_data(JH15_ASV) <- design 
sample_data(JH15_ASV)

##################################################################################
##Constrained ordination
####################################################################################

##We are going to subset the phyloseq object 
JH15_genotype_HID2<- subset_samples(JH15_ASV, Genotype == "Barke" |Genotype == "52"|Genotype == "17") 
JH15_genotype_HID2

#Constrained ordination: constrained for Description
#rarefied data
JH15_CAP_genotype_HID2 <- ordinate(JH15_genotype_HID2, "CAP", "bray", ~ Description)
plot_ordination(JH15_genotype_HID2, JH15_CAP_genotype_HID2, color = "Description")

#Assign shapes to Soil type and color to Sample type
p=plot_ordination(JH15_genotype_HID2, JH15_CAP_genotype_HID2 , color = "Genotype", shape= "Genotype")
p = p + geom_point(size = 6, alpha = 0.80)
p = p + scale_colour_manual(values = c("#56B4E9","#E69F00","#0072B2"))#"#56B4E9"
p + ggtitle("CAP 16S data, Bray distance-rarefied samples")

#Extract the mapping file
design_rhizosphere_HID2 <- design[colnames(otu_table(JH15_genotype_HID2)), ]

#BC distance
BC <- phyloseq::distance(JH15_genotype_HID2, "bray")

adonis(BC ~ Genotype, data= design_rhizosphere_HID2, permutations = 5000)

adonis(BC ~ Sampling_date*Genotype , data=design_rhizosphere_HID2, permutation =5000)

######################################################################################
#DEseq 
######################################################################################
##Subset the samples
JH15_17_52 <- subset_samples(JH15_ASV, Genotype == "17"|Genotype == "52"| Genotype == "Barke"|Genotype == "Bulk") 
JH15_17_52

#extract count data : 
JH15_counts_17_52 <- otu_table(JH15_17_52)
countData_17_52 = as.data.frame(JH15_counts_17_52)
colnames(JH15_counts_17_52)

#the design file containing sample information
colData_17_52 = design[colnames(JH15_counts_17_52), ]
rownames(colData_17_52)
colData_17_52

#construct a DESeq dataset combining count data and sample information
JH15_cds_17_52 <- DESeqDataSetFromMatrix(countData =countData_17_52, colData=colData_17_52 , design= ~Description)
JH15_cds_17_52

#execute the differential count analysis with the function DESeq 
JH15_cds_test_17_52 <- DESeq(JH15_cds_17_52, fitType="local", betaPrior = FALSE) 
JH15_cds_test_17_52

#saveRDS(JH15_cds_test_17_52, file = "JH15_cds_test_17_52.rds")
JH15_cds_test_17_52<-readRDS("JH15_cds_test_17_52.rds")

#Disclaimer to identify OTU enriched in B124.52 vs B.124.17
########################################################################
#Define the OTUs significantly enriched in the rhizosphere samples Bulk vs 124.52 (_17 to differentialte preovious comparisons)

rhizo_B124.52_a <- results(JH15_cds_test_17_52, contrast = c("Description","BULK","B124.HID2")) 

#Inspect a result file
rhizo_B124.52_a

#Extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_B124.52_a_FDR005 <- rhizo_B124.52_a[(rownames(rhizo_B124.52_a)[which(rhizo_B124.52_a$padj <0.05)]), ]

#Enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_B124.52_a_enriched <-  rhizo_B124.52_a[(rownames(rhizo_B124.52_a)[which(rhizo_B124.52_a$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_B124.52_a_enriched_FDR005 <- intersect(rownames(rhizo_B124.52_a_FDR005), rownames(rhizo_B124.52_a_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_B124.52_a_enriched_FDR005)


#Identify the ones differential enriched between Barke and Bulk soil
###########################################################################################
#Define the OTUs significantly enriched in the rhizosphere samples
rhizo_Barke_a <- results(JH15_cds_test_17_52, contrast = c("Description", "BARKE", "BULK"))

#Inspect a result file
rhizo_Barke_a

#Extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_Barke_a_FDR005 <- rhizo_Barke_a[(rownames(rhizo_Barke_a)[which(rhizo_Barke_a$padj <0.05)]), ]

#Enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_Barke_a_enriched <-  rhizo_Barke_a[(rownames(rhizo_Barke_a)[which(rhizo_Barke_a$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_Barke_a_enriched_FDR005 <- intersect(rownames(rhizo_Barke_a_FDR005), rownames(rhizo_Barke_a_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_Barke_a_enriched_FDR005)

#Identify the ones differential enriched between B.124.17 and Bulk soil
################################################################################
#Define the OTUs significantly enriched in the rhizosphere samples
rhizo_B124.Barke_a <- results(JH15_cds_test_17_52, contrast = c("Description", "B124.BARKE", "BULK"))

#Inspect a result file
rhizo_B124.Barke_a

#Extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_B124.Barke_a_FDR005 <- rhizo_B124.Barke_a[(rownames(rhizo_B124.Barke_a)[which(rhizo_B124.Barke_a$padj <0.05)]), ]

#Enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_B124.Barke_a_enriched <-  rhizo_B124.Barke_a[(rownames(rhizo_B124.Barke_a)[which(rhizo_B124.Barke_a$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_B124.Barke_a_enriched_FDR005 <- intersect(rownames(rhizo_B124.Barke_a_FDR005), rownames(rhizo_B124.Barke_a_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_B124.Barke_a_enriched_FDR005)

#Identify the ones differential enriched between Barke and B124.HID2_52
###############################################################################
#Define the OTUs significantly enriched in the rhizosphere samples
B124.52_Barke_a <- results(JH15_cds_test_17_52, contrast = c("Description", "BARKE", "B124.HID2"))

#inspect a result file
B124.52_Barke_a

#Extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
B124.52_Barke_a_FDR005 <- B124.52_Barke_a[(rownames(B124.52_Barke_a)[which(B124.52_Barke_a$padj <0.05)]), ]
length(B124.52_Barke_a_FDR005)
#Enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
B124.52_Barke_a_enriched <-  B124.52_Barke_a[(rownames(B124.52_Barke_a)[which(B124.52_Barke_a$log2FoldChange > 0)]), ]
length (B124.52_Barke_a_enriched)
#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
B124.52_a_enriched_FDR005_vs_Barke <- intersect(rownames(B124.52_Barke_a_FDR005), rownames(B124.52_Barke_a_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(B124.52_a_enriched_FDR005_vs_Barke)

JH15_no_plants_dada2[B124.52_a_enriched_FDR005_vs_Barke, ]

#Are these OTUs enriched from soil?
B124.52_a_enriched_FDR005_vs_Barke_rhizo <- intersect(rhizo_B124.52_a_enriched_FDR005, B124.52_a_enriched_FDR005_vs_Barke)
length(B124.52_a_enriched_FDR005_vs_Barke_rhizo)
B124.52_a_enriched_FDR005_vs_Barke_rhizo
JH15_no_plants_dada2[B124.52_a_enriched_FDR005_vs_Barke_rhizo, ]

####End
##############################################################################################################################################

#######################################################################
#######################################################################
#Supplementary_figure 7 Alpha diversity calculations 
#######################################################################
#######################################################################

JH15_alpha_rare_52_17 <-  estimate_richness(JH15_17_52, measures = c("Observed", "Shannon", "Chao1")) 
JH15_alpha_rare_52_17 

JH15_17_52_otu_table<-otu_table(JH15_17_52)
#Data frame Genotype_Description 
colData_17_52 = design[colnames(JH15_17_52_otu_table), ]
rownames(colData_17_52)
colData_17_52


#Description 
design_genotype  <- as.data.frame(colData_17_52[, 1]) 
rownames(design_genotype) <- rownames(colData_17_52) 
colnames(design_genotype) <- c("Genotype") 
design_genotype  

###############################
#### Alpha diversity OBSERVED 
##############################
#Observed ASVs 
JH15_alpha_rare_52_17_Observed <- as.data.frame(JH15_alpha_rare_52_17[ ,1]) 
rownames(JH15_alpha_rare_52_17_Observed) <- rownames(JH15_alpha_rare_52_17) 
colnames(JH15_alpha_rare_52_17_Observed) <- c("Observed") 

#Combine the dataset sample description and Observed OTUs 
JH15_alpha_rare_52_17_Observed_TD <- cbind(design_genotype , JH15_alpha_rare_52_17_Observed) 
JH15_alpha_rare_52_17_Observed_TD <- as.data.frame(JH15_alpha_rare_52_17_Observed_TD)
JH15_alpha_rare_52_17_Observed_TD$Genotype


#Order the levels according to a defined order 
JH15_alpha_rare_52_17_Observed_TD$Genotype <- ordered(JH15_alpha_rare_52_17_Observed_TD$Genotype, levels=c("Bulk", "Barke", "17","52"))  

#Plotting
with(JH15_alpha_rare_52_17_Observed_TD, boxplot(Observed~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", main = "ASVs observed richness", col=c("grey","#0072B2","#56B4E9","#E69F00", "#D55E00")))
with(JH15_alpha_rare_52_17_Observed_TD, stripchart(Observed~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Observed_OTUs_stats <- aov(Observed ~ Genotype, data = JH15_alpha_rare_52_17_Observed_TD) 
summary(Observed_OTUs_stats) 

#Tukey test 
TukeyHSD(Observed_OTUs_stats, ordered = TRUE) 

###########################
#### Alpha diversity CHAO1 
##########################
#Chao1 ASVs 
JH15_alpha_rare_52_17_Chao1 <- as.data.frame(JH15_alpha_rare_52_17[ ,2]) 
rownames(JH15_alpha_rare_52_17_Chao1) <- rownames(JH15_alpha_rare_52_17) 
colnames(JH15_alpha_rare_52_17_Chao1) <- c("Chao1") 


#Combine the dataset sample description and Chao1 OTUs 
JH15_alpha_rare_52_17_Chao1_TD <- cbind(design_genotype , JH15_alpha_rare_52_17_Chao1) 
JH15_alpha_rare_52_17_Chao1_TD <- as.data.frame(JH15_alpha_rare_52_17_Chao1_TD) 
JH15_alpha_rare_52_17_Chao1_TD$Genotype


#Order the levels according to a defined order 
JH15_alpha_rare_52_17_Chao1_TD$Genotype <- ordered(JH15_alpha_rare_52_17_Chao1_TD$Genotype, levels=c("Bulk", "Barke", "17","52"))  

#Plotting
with(JH15_alpha_rare_52_17_Chao1_TD, boxplot(Chao1~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", main = "ASVs Chao richness", col=c("grey","#0072B2","#56B4E9","#E69F00", "#D55E00")))
with(JH15_alpha_rare_52_17_Chao1_TD, stripchart(Chao1~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))

#ANOVA 
Chao1_OTUs_stats <- aov(Chao1 ~ Genotype, data = JH15_alpha_rare_52_17_Chao1_TD) 
summary(Chao1_OTUs_stats) 

#Tukey test 
TukeyHSD(Chao1_OTUs_stats, ordered = TRUE) 

#############################
#### Alpha diversity SHANNON
#############################
#Shannon ASvs 
JH15_alpha_rare_52_17_Shannon <- as.data.frame(JH15_alpha_rare_52_17[ ,4]) 
rownames(JH15_alpha_rare_52_17_Shannon) <- rownames(JH15_alpha_rare_52_17) 
colnames(JH15_alpha_rare_52_17_Shannon) <- c("Shannon") 


#Combine the dataset sample description and Shannon OTUs 
JH15_alpha_rare_52_17_Shannon_TD <- cbind(design_genotype , JH15_alpha_rare_52_17_Shannon) 
JH15_alpha_rare_52_17_Shannon_TD <- as.data.frame(JH15_alpha_rare_52_17_Shannon_TD) 
JH15_alpha_rare_52_17_Shannon_TD$Genotype

#Order the levels according to a defined order 
JH15_alpha_rare_52_17_Shannon_TD$Genotype <- ordered(JH15_alpha_rare_52_17_Shannon_TD$Genotype, levels=c("Bulk", "Barke", "17","52"))  

#Plotting
with(JH15_alpha_rare_52_17_Shannon_TD, boxplot(Shannon~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", main = "ASVs Shannon richness",col=c("grey","#0072B2","#56B4E9","#E69F00", "#D55E00")))
with(JH15_alpha_rare_52_17_Shannon_TD, stripchart(Shannon~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Shannon_OTUs_stats <- aov(Shannon ~ Genotype, data = JH15_alpha_rare_52_17_Shannon_TD) 
summary(Shannon_OTUs_stats) 

#Tukey test 
TukeyHSD(Shannon_OTUs_stats, ordered = TRUE)

############################################################################################################################################
##End


############################################################
############################################################
###Supplementary figure 8: 16S and ITS DNA copies
############################################################
#########################################################

##The CT values of the qPCR were calculated as follow: 

#1) Linear extrapolation of the samples CT values in the standard curve (CT standards and Ln of concentrations (fg)), with the formula: y=ax+b, where y=CT value of the sample
#Y(x) = Y(1)+ (x- x(1)/x(2)-x(1)) * (Y(2) - Y(1)) in Excel

#2) Transformed in copy numbers:
    ## Bacteria 16S= Genome copy # = DNA (g) / (g_to_bp const. x genome size) With DNA (g) = 20 ng g to bp const = 1.096 x 10^-21 g genome size = 4.6 x 10^6 bp 20 x 10^-9 g / (1.096 x 10^-21 g/bp x 4.6 x 10^6 bp) = 4.0 x 10^6 copies in 20 ng
    ##Fungi ITS=  Genome copy # = DNA (g) / (g_to_bp const. x genome size) With DNA (g) = 20 ng g to bp const = ITS copy#=20ng* 1.079 x 10^-12 genome size#=12.1*10^6  20 ng / ( 1.079 x 10^-12 ng/bp x 12.1 x 10^6 bp) = 1531874 copies in 20 ng

#3) Log transformed

qJH23 <-(read.delim("JH15_JH23_quantification_copy_numbers.txt", sep = "\t"))

###########################################
###############Plotting Bacteria 16S copies
###########################################

#order the factor
qJH23$Sample <-  ordered(qJH23$Sample, levels=c( "Bulk","Barke", "Intr_17", "Intr_52"))

p <-ggplot(qJH23, aes(x=Sample, y=Log_16S_copies, fill=Sample)) + geom_boxplot() + ggtitle("Bacterial Log 16S copy numbers")+ylim(5,25)
p + geom_jitter( size=5,shape=21, position=position_jitter(0.01))+ scale_fill_manual(values = c("#000000","#0072B2","#56B4E9", "#E69F00", "#D55E00","#D55E00"))

##################################################################

shapiro.test(qJH23$Log_16S_copies)

#Stats Bacterial concentrations
kruskal.test(Log_16S_copies~Sample , data=qJH23)

posthoc.kruskal.dunn.test (x= qJH23$Log_16S_copies, g= qJH23$Sample, p.adjust.method="BH")

###################################################################
###############plotting concentrations Log Fungal concentrations
##################################################################

qJH23$Sample <-  ordered(qJH23$Sample, levels=c( "Bulk","Barke", "Intr_17", "Intr_52"))


q <-ggplot(qJH23, aes(x=Sample, y=Log_ITS_copies, fill=Sample)) + geom_boxplot() + ggtitle(" Log Fungal ITS copies")+ylim(5,25)
q + geom_jitter( size=5,shape=21, position=position_jitter(0))+ scale_fill_manual(values = c ("#000000","#0072B2","#56B4E9", "#E69F00", "#D55E00","#D55E00"))

##################################################################
#Stats Fungal ratios
shapiro.test(qJH23$Log_ITS_copies)

kruskal.test(Log_ITS_copies~Sample , data=qJH23)

posthoc.kruskal.dunn.test (x= qJH23$Log_ITS_copies, g= qJH23$Sample, p.adjust.method="BH")


##End
################################################################################################################################################



















