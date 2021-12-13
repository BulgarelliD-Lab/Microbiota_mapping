#####################################################################################
#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:
#  Figure 2
#  Revision 11/21 
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
library("qtl2")
library("devtools")
#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

##################################################################################
## Figure 2
##################################################################################

#import the file: 
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


########################################################################################
JH07_parentals <- subset_samples(JH07_ASV, Description != "BC1S3")

JH07_parentals_rhizo <- subset_samples(JH07_ASV, Description == "Barke"| Description == "HID144")

#calculation of differentially abundant taxa between the parental lines using Deseq2
######################################################################################

#extract count data 
JH07_counts <- otu_table(JH07_parentals)
countData = as.data.frame(JH07_counts)
colnames(JH07_counts)

#the design file containing sample information
colData = design[colnames(JH07_counts), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH07_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH07_cds_test <- DESeq(JH07_cds, fitType="local", betaPrior = FALSE)

#saveRDS(JH07_cds_test, file="JH07_cds_test_ASVs")
JH07_cds_test<-readRDS(file="JH07_cds_test_ASVs")
levels(JH07_cds_test$Description)


#define the OTUs significantly enriched in the rhizosphere samples
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

##########################################################################################
#Mapping exercise
###########################################################################################
#extract the count data for the whole population
#extract count data 
#remove soil samples
JH07_population <- subset_samples(JH07_ASV, Description != "Bulk")
JH07_counts_population <- otu_table(JH07_population)
JH07_counts_population = as.data.frame(JH07_counts_population)
colnames(JH07_counts_population)
#create a new design file
JH07_design_population <- design[colnames(JH07_counts_population), ]
#replace the names in the OTU table with genotype IDs
colnames(JH07_counts_population) <- JH07_design_population$Genotype
JH07_counts_population[1:10, ]
#generate average values for the OTU counts
JH07_counts_population_average <- as.data.frame(do.call(cbind, by(t(JH07_counts_population),INDICES=names(JH07_counts_population),FUN=colMeans)))
JH07_counts_population_average[1:5, ]

#extract the count information for the Elite enriched versus wild
JH07_counts_population_average_HID_enriched <- JH07_counts_population_average[rownames(HID144_enriched_Barke), ]
JH07_counts_population_average_Barke_enriched <- JH07_counts_population_average[rownames(Barke_enriched_HID144), ]

dim(JH07_counts_population_average_HID_enriched)
dim(JH07_counts_population_average_Barke_enriched)

#create the dataset for the mapping exercise

#retrieve taxonomy information
JH07_HID144_ASVs_enriched<-JH07_no_plants_dada2[rownames(JH07_counts_population_average_HID_enriched), ]
JH07_Barke_ASVs_enriched<-JH07_no_plants_dada2[rownames(JH07_counts_population_average_Barke_enriched), ]

JH07_HID144_ASVs_ASV_enriched<- cbind(JH07_counts_population_average_HID_enriched,JH07_HID144_ASVs_enriched)
JH07_Barke_ASVs_ASV_enriched<- cbind(JH07_counts_population_average_Barke_enriched,JH07_Barke_ASVs_enriched)
#write.table(JH07_HID144_ASVs_ASV_enriched, file="JH07_merged_HID144_Silva138_ASV_enriched_Mapping_wo_N1", sep="\t")
#write.table(JH07_Barke_ASVs_ASV_enriched, file="JH07_merged_Barke_Silva138_ASV_enriched_Mapping_wo_N1", sep="\t")

#####################################################################################################################
##QTL mapping ASV level
#####################################################################################################################

###Read the file and define the population type
library(qtl)
library(devtools)

##Read.cross objects:https://www.rdocumentation.org/packages/qtl/versions/1.46-2/topics/read.cross

wild <- read.cross("csv", ".", "JH07_merged_HID144-Barke_Silva138_ASV_enriched_Mapping_0621.csv", BC.gen = 1, F.gen = 3)
###As the map position comes from genetic map wherein several markers map to same location, we might get error message but that's fine.
wild <- jittermap(wild)
summary(wild)
plot.map(wild)

###calculate the genotypic probability with almost default parameters
wild1<-calc.genoprob(wild, step=2.5,error.prob=0.01)
###similarity of the genotype
wild2<-sim.geno(wild, step=2.5, n.draws=64, error.prob=0.01)
###Different methods of analysis see the manual for details, specify the phenotype column
#Change here the number of phenotypes
out.em1 <- scanone(wild1, pheno.col =1:36) #number of phenotypes or taxa = 36
###get the summary of the analysis
summary(out.em1)

###Another method
out.hk1 <- scanone(wild1, pheno.col =1:36, method="hk")
###get the summary of the analysis
summary(out.hk1)
###Another method
out.imp1 <- scanone(wild2, pheno.col = 1:36, method="imp")
###get the summary of the analysis
summary(out.imp1)


###writing the results in a csv file
#####################################
#write.csv(out.em1, "JH07_merged_ASVs_HID144_Barke_Silva138_mapping_0621_out_em1_wo_het.csv")
#write.csv(out.hk1, "JH07_merged_ASVs_HID144_Silva138_mapping_0321_out_hk1.csv")
#write.csv(out.imp1, "JH07_merged_ASVs_HID144_Silva138_mapping_0321_out_imp1.csv")
#output LOD thresholds 


###################################################################################
###2. to get a LOD genome-wide significance threshold or genome-scan-adjusted p-values
####################################################################################
perm <- scanone(wild1, method="em", n.perm=1000, pheno.col =1:36)
plot(perm)

#Significance thresholds may be obtained via the summary() function:
summary(perm)
summary(perm, alpha=c(0.05, 0.2))

#Having significance thresholds and p-values calculated automatically:

summary(out.em1, perms=perm, alpha=0.2, pvalues=TRUE)

#Significance thresholds were calculated by permutation, with thresholds of 20% and 5% genome wide significance 

#################################
##3. Test for epistasis
##################################
summary.scantwo(out22.em)

out22.em <-scantwo(wild1)

plot(out22.em)

summary(out22.em)


operm2 <-scantwo(wild1,method="hk",n.perm=100, pheno.col =1:36)

summary(operm2, alpha=c(0.05,0.20))


summary(out22.em)#NO evidence for epistasis

plot(out2.em, chr=c(3,5))

max(out22.em)


###################################
##4. Fitting the model
###################################
#You can use fitqtl. It's marked as "%var" in the summary. The %var assigned to individual QTL are not trustworthy in the case of linked QTL. 

##R2 calculation per taxa (phenotype)
#########################################

###Enchiched in Elite
#Lysobacter
qtl_Lysobacter <- makeqtl(wild2, chr=c("2H","2H","5H","6H","7H"), pos=c(37.8,126.8,89.5,89.35,22.7))
summary(out.qtl_Lysobacter <- fitqtl(wild2, pheno.col=30, qtl=qtl_Lysobacter , formula=y~Q1+Q2+Q3+Q4+Q5,method="imp", get.ests=TRUE))
##No significant QTLs

##Streptomyces4
qtl_Streptomyces4 <- makeqtl(wild2, chr=c("2H","4H","7H"), pos=c(37.8,72.5,75.85))
summary(out.qtl_Streptomyces4 <- fitqtl(wild2, pheno.col=36, qtl=qtl_Streptomyces4 , formula=y~Q1+Q2+Q3,method="imp", get.ests=TRUE))

##Drop no significant QTLs
qtl_Streptomyces4 <- makeqtl(wild2, chr=c("4H"), pos=c(72.5))
summary(out.qtl_Streptomyces4 <- fitqtl(wild2, pheno.col=36, qtl=qtl_Streptomyces4, formula=y~Q1,method="imp", get.ests=TRUE))
#R2 37.80721 p=1.803104e-05
qtl_Streptomyces4 <- makeqtl(wild2, chr=c("2H"), pos=c(37.8))
summary(out.qtl_Streptomyces4 <- fitqtl(wild2, pheno.col=36, qtl=qtl_Streptomyces4, formula=y~Q1,method="imp", get.ests=TRUE))

qtl_Streptomyces4 <- makeqtl(wild2, chr=c("7H"), pos=c(75.85))
summary(out.qtl_Streptomyces4 <- fitqtl(wild2, pheno.col=36, qtl=qtl_Streptomyces4, formula=y~Q1,method="imp", get.ests=TRUE))
##No significant QTLs

##Norcardia
qtl_Norcardia <- makeqtl(wild2, chr=c("7H"), pos=c(140.6))
summary(out.qtl_Norcardia <- fitqtl(wild2, pheno.col=31, qtl=qtl_Norcardia, formula=y~Q1,method="imp", get.ests=TRUE))
##No significant QTLs


############################################
###Enchiched in Wild
###########################################
##Pedobacter2
qtl_Pedobacter2 <- makeqtl(wild2, chr=c("2H","6H"), pos=c(92.95,0))
summary(out.qtl_Pedobacter2 <- fitqtl(wild2, pheno.col=25, qtl=qtl_Pedobacter2 , formula=y~Q1+Q2,method="imp", get.ests=TRUE))
##Drop no significant QTLs
qtl_Pedobacter2 <- makeqtl(wild2, chr=c("2H"), pos=c(92.95))
summary(out.qtl_Pedobacter2 <- fitqtl(wild2, pheno.col=25, qtl=qtl_Pedobacter2 , formula=y~Q1,method="imp", get.ests=TRUE))

qtl_Pedobacter2 <- makeqtl(wild2, chr=c("6H"), pos=c(0))
summary(out.qtl_Pedobacter2 <- fitqtl(wild2, pheno.col=25, qtl=qtl_Pedobacter2 , formula=y~Q1,method="imp", get.ests=TRUE))

#Vicinamibacterales #from 9.6 to 12.5/
qtl_Vicinamibacterales <- makeqtl(wild2, chr=c("3H","5H"), pos=c(9.6,96.6))
summary(out.qtl_Vicinamibacterales <- fitqtl(wild2, pheno.col=16, qtl=qtl_Vicinamibacterales , formula=y~Q1+Q2,method="imp", get.ests=TRUE))

##Sorangium
qtl_sorangium <- makeqtl(wild2, chr=c("3H"), pos=c(38.75)) #based on imputed genotypes
summary(out.qtl_sorangium<- fitqtl(wild2, pheno.col=18, qtl=qtl_sorangium , formula=y~Q1,method="hk", get.ests=TRUE))

##Holophaga
qtl_holophaga <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.qtl_holophaga <- fitqtl(wild2, pheno.col=15, qtl=qtl_holophaga , formula=y~Q1,method="imp", get.ests=TRUE))

##Tahibacter
qtl_tahibacter <- makeqtl(wild2, chr=c("3H","3H"), pos=c(38.75,49.9))
summary(out.qtl_tahibacter <- fitqtl(wild2, pheno.col=21, qtl=qtl_tahibacter , formula=y~Q1+Q2,method="imp", get.ests=TRUE))

##Drop no significant QTLs
qtl_tahibacter <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.qtl_tahibacter <- fitqtl(wild2, pheno.col=21, qtl=qtl_tahibacter, formula=y~Q1,method="imp", get.ests=TRUE))

##Variovorax
qtl_Variovorax <- makeqtl(wild2, chr=c("3H","3H"), pos=c(38.75,47.3))
summary(out.qtl_Variovorax <- fitqtl(wild2, pheno.col=11, qtl=qtl_Variovorax , formula=y~Q1+Q2,method="imp", get.ests=TRUE))

##Drop no significant QTLs
qtl_Variovorax <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.qtl_Variovorax <- fitqtl(wild2, pheno.col=11, qtl=qtl_Variovorax , formula=y~Q1,method="imp", get.ests=TRUE))

#R2 63.33367 p= 9.509349e-11


##Bdellovibrio
qtl_Bdellovibrio <- makeqtl(wild2, chr=c("4H"), pos=c(91.7))
summary(out.qtl_Bdellovibrio <- fitqtl(wild2, pheno.col=13, qtl=qtl_Bdellovibrio , formula=y~Q1,method="imp", get.ests=TRUE))

##Stenotrophomonas 
qtl_Stenotrophomonas <- makeqtl(wild2, chr=c("4H", "6H", "7H"), pos=c(94.1,0, 128.3))
summary(out.qtl_Stenotrophomonas <- fitqtl(wild2, pheno.col=7, qtl=qtl_Stenotrophomonas , formula=y~Q1+Q2+Q3,method="imp", get.ests=TRUE))

#drop No significant

qtl_Stenotrophomonas <- makeqtl(wild2, chr=c("4H"), pos=c(94.1))
summary(out.qtl_Stenotrophomonas <- fitqtl(wild2, pheno.col=7, qtl=qtl_Stenotrophomonas , formula=y~Q1,method="imp", get.ests=TRUE))


########################################################
##Confidence intervals for the significant QTLs
#######################################################
#To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL of interest, type:


##Vicinamibacterales
bayesint(out.em1, chr="2H", out.2qtl, prob=0.95, lodcolumn=16, expandtomarkers=FALSE)
bayesint(out.em1, chr="5H", out.2qtl, prob=0.95, lodcolumn=16, expandtomarkers=FALSE)

###Variovorax
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=11, expandtomarkers=FALSE)
###Sorangium
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=18, expandtomarkers=FALSE)

###Holophaga
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=15, expandtomarkers=FALSE)

##Tahibacter
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=21, expandtomarkers=FALSE)

##Stenothrophomonas
bayesint(out.em1, chr="4H", out.2qtl, prob=0.95, lodcolumn=7, expandtomarkers=FALSE)

##Streptomyces4
bayesint(out.em1, chr="4H", out.2qtl, prob=0.95, lodcolumn=36, expandtomarkers=FALSE)


##################################################
##Marker regresion
###################################################
##Allele composition at certain marker location:Plot the phenotype values versus the genotypes at a marker or markers. 

##Streptomyces4
Streptomyces4plotPXG<-plotPXG(wild, "BOPA1_8653_475", pheno.col=36, jitter=1, infer=T)
Streptomyces4aov<-aov(Streptomyces4plotPXG$pheno~Streptomyces4plotPXG$BOPA1_8653_475) 
summary(Streptomyces4aov)
TukeyHSD(Streptomyces4aov)

##Vicinamibacterales

VicinamibacteralesplotPXG<-plotPXG(wild, "SCRI_RS_108971", pheno.col=16, jitter=1, infer=T)
shapiro.test(VicinamibacteralesplotPXG$pheno)
Vicinamibacteraleskruskal<-kruskal.test(VicinamibacteralesplotPXG$pheno~VicinamibacteralesplotPXG$SCRI_RS_108971) 
Vicinamibacteraleskruskal

dunn.test(x= VicinamibacteralesplotPXG$pheno, g= VicinamibacteralesplotPXG$SCRI_RS_108971)


VicinamibacteralesplotPXG<-plotPXG(wild, "SCRI_RS_206565", pheno.col=16, jitter=1, infer=T)

Vicinamibacteralesaov<-aov(VicinamibacteralesplotPXG$pheno~VicinamibacteralesplotPXG$SCRI_RS_206565) 
summary(Vicinamibacteralesaov)
TukeyHSD(Vicinamibacteralesaov)
##Sorangium
SorangiumplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=18, jitter=1, infer=T)

shapiro.test(SorangiumplotPXG$pheno)
Sorangiumkruskal<-kruskal.test(SorangiumplotPXG$pheno~SorangiumplotPXG$SCRI_RS_154747) 
Sorangiumkruskal
Sorangiumaov<-aov(SorangiumplotPXG$pheno~SorangiumplotPXG$SCRI_RS_154747) 
summary(Sorangiumaov)
TukeyHSD(Sorangiumaov)
#Varivorax

VarivoraxplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=11, jitter=1, infer=T)

Varivoraxaov<-aov(VarivoraxplotPXG$pheno~VarivoraxplotPXG$SCRI_RS_154747) 
summary(Varivoraxaov)
TukeyHSD(Varivoraxaov)

##Holophaga
HolophagaplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=15, jitter=1, infer=T)#,

Holophagaaov<-aov(HolophagaplotPXG$pheno~HolophagaplotPXG$SCRI_RS_154747) 
summary(Holophagaaov)
TukeyHSD(Holophagaaov)

##Tahibacter
TahibacterplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=21, jitter=1, infer=T)#,

Tahibacteraov<-aov(TahibacterplotPXG$pheno~TahibacterplotPXG$SCRI_RS_154747) 
summary(Tahibacteraov)
TukeyHSD(Tahibacteraov)
###Stenotrophomonas
StenotrophomonasplotPXG<-plotPXG(wild, "SCRI_RS_188944", pheno.col=7, jitter=1, infer=T)#,

Stenotrophomonasaov<-aov(StenotrophomonasplotPXG$pheno~StenotrophomonasplotPXG$SCRI_RS_188944) 
summary(Stenotrophomonasaov)
TukeyHSD(Stenotrophomonasaov)


############################################################################################################################
##End

#####################################################################################################################
##QTL mapping Genus level
#####################################################################################################################


###Read the file and define the population type
library(qtl)
library(devtools)

##Read.cross objects:https://www.rdocumentation.org/packages/qtl/versions/1.46-2/topics/read.cross

wild <- read.cross("csv", ".", "JH07_merged_HID144-Barke_Silva138_Genus_enriched_Mapping.csv", BC.gen = 1, F.gen = 3)
###As the map position comes from genetic map wherein several markers map to same location, we might get error message but that's fine.
wild <- jittermap(wild)
summary(wild)
plot.map(wild)


###calculate the genotypic probability with almost default parameters
wild1<-calc.genoprob(wild, step=2.5,error.prob=0.01)
###similarity of the genotype
wild2<-sim.geno(wild, step=2.5, n.draws=64, error.prob=0.01)
###Different methods of analysis see the manual for details, specify the phenotype column
#Change here the number of phenotypes
out.em1 <- scanone(wild1, pheno.col =1:18) #number of phenotypes or taxa = 18
###get the summary of the analysis
summary(out.em1)#I think here is giving feature only for the first  column taxa Pedobacter
max(out.em1)


###writing the results in a csv file, here just an example
#write.csv(out.em1, "JH07_merged_Genus_Barke_Silva138_mapping_0321_out_em1.csv")
###Another method
out.hk1 <- scanone(wild1, pheno.col =1:18, method="hk")
###get the summary of the analysis
summary(out.hk1)
###Another method
out.imp1 <- scanone(wild2, pheno.col = 1:18, method="imp")
###get the summary of the analysis
summary(out.imp1)

#write.csv(out.hk1, "JH07_merged_Genus_Barke_Silva138_mapping_0321_out_em1.csv")
#write.csv(out.imp1, "JH07_merged_Genus_Barke_Silva138_mapping_0321_out_em1.csv")
#output LOD thresholds 

###################################################################################
###2. to get a LOD genome-wide significance threshold or genome-scan-adjusted p-values
####################################################################################

###to get a LOD genome-wide significance threshold or genome-scan-adjusted p-values
perm <- scanone(wild1, method="em", n.perm=1000, pheno.col =1:18)
plot(perm)


#Significance thresholds may be obtained via the summary() function:
summary(perm)
summary(perm, alpha=c(0.05, 0.2))

#Having significance thresholds and p-values calculated automatically:
summary(out.em1, perms=perm, alpha=0.2, pvalues=TRUE)


#################################
##3. Test for epistasis
##################################

out22.em <-scantwo(wild1, pheno.col =1:18)

plot(out22.em)

summary(out22.em)


operm2 <-scantwo(wild1,method="hk",n.perm=100)

summary(operm2, alpha=c(0.05,0.20))


summary(out22.em)#NO evidence for epistasis

plot(out2.em, chr=c(3,5))

max(out22.em)

########################################################
##Confidence intervals
#######################################################
#To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL on chromosome 3, type:
#Ramlibacter

bayesint(out.em1, chr="2H", out.2qtl, prob=0.95, lodcolumn=3, expandtomarkers=FALSE)
##Variovorax
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=7, expandtomarkers=FALSE)

###Sorangium
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=13, expandtomarkers=FALSE)

##Rhodanobacter
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=4, expandtomarkers=FALSE)

###Holophaga
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=15, expandtomarkers=FALSE)

##Tahibacter
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=21, expandtomarkers=FALSE)

##Microbacterium
bayesint(out.em1, chr="5H", out.2qtl, prob=0.95, lodcolumn=10, expandtomarkers=FALSE)

##Streptomyces
bayesint(out.em1, chr="7H", out.2qtl, prob=0.95, lodcolumn=18, expandtomarkers=FALSE)

###################################
##4. Fitting the model
###################################

#You can use fitqtl. It's marked as "%var" in the summary. The %var assigned to individual QTL are not trustworthy in the case of linked QTL. 


##R2 calculation per taxa (phenotype)
##########################################
##Enriched in wild
############################
#Ramlibacter
qtl_variovorax <- makeqtl(wild2, chr=c("2H"), pos=c(123.5)) #based on imputed genotypes

summary(out.2qtl <- fitqtl(wild2, pheno.col=3, qtl=qtl_variovorax, formula=y~Q1,method="imp", get.ests=TRUE))

##Variovorax
qtl_variovorax <- makeqtl(wild2, chr=c("3H"), pos=c(38.75)) #based on imputed genotypes
summary(out.2qtl <- fitqtl(wild2, pheno.col=7, qtl=qtl_variovorax, formula=y~Q1,method="imp", get.ests=TRUE))


##Sorangium
qtl_sorangium <- makeqtl(wild2, chr=c("1H","3H"), pos=c(97.5,38.75)) #based on imputed genotypes

summary(out.2qtl <- fitqtl(wild2, pheno.col=13, qtl=qtl_sorangium , formula=y~Q1+Q2,method="imp", get.ests=TRUE))
##Drop no significant QTLs
qtl_sorangium <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.2qtl <- fitqtl(wild2, pheno.col=13, qtl=qtl_sorangium , formula=y~Q1,method="imp", get.ests=TRUE))

##Rhodanobacter
qtl_rhodanobacter <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.2qtl <- fitqtl(wild2, pheno.col=4, qtl=qtl_rhodanobacter , formula=y~Q1,method="imp", get.ests=TRUE))

##Holophaga
qtl_holophaga <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.2qtl <- fitqtl(wild2, pheno.col=12, qtl=qtl_holophaga , formula=y~Q1,method="imp", get.ests=TRUE))

##Tahibacter
qtl_tahibacter <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.2qtl <- fitqtl(wild2, pheno.col=15, qtl=qtl_tahibacter , formula=y~Q1,method="imp", get.ests=TRUE))

##P3OB.42
qtl_P3OB.42 <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.2qtl <- fitqtl(wild2, pheno.col=11, qtl=qtl_P3OB.42 , formula=y~Q1,method="imp", get.ests=TRUE))

##Enriched in elite
######################
##Microbacterium
qtl_Microbacterium <- makeqtl(wild2, chr=c("5H"), pos=c(118.1))
summary(out.2qtl <- fitqtl(wild2, pheno.col=10, qtl=qtl_Microbacterium  , formula=y~Q1,method="imp", get.ests=TRUE))

####Streptomyces
qtl_Streptomyces <- makeqtl(wild2, chr=c("7H"), pos=c(133.9))
summary(out.2qtl <- fitqtl(wild2, pheno.col=18, qtl=qtl_Streptomyces  , formula=y~Q1,method="imp", get.ests=TRUE))


##################################################
##Marker regresion
###################################################
##Allele composition at certain marker location:Plot the phenotype values versus the genotypes at a marker or markers. 
##T-test is used in BC populations
##Ramlibacter

RamlibacterplotPXG<-plotPXG(wild, "BOPA1_285_2932", pheno.col=3, jitter=1, infer=T)

Ramlibacteraov<-aov(RamlibacterplotPXG$pheno~RamlibacterplotPXG$BOPA1_285_2932) 
summary(Ramlibacteraov)
TukeyHSD(Ramlibacteraov) 
############

#Variovorax
VariovoraxplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=7, jitter=1, infer=T)
Variovoraxaov<-aov(VariovoraxplotPXG$pheno~VariovoraxplotPXG$SCRI_RS_154747) 
summary(Variovoraxaov) 
#t.test(VariovoraxplotPXG$pheno~VariovoraxplotPXG$SCRI_RS_154747) 
TukeyHSD(Variovoraxaov)

##Sorangium
SorangiumplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=13, jitter=1, infer=T)
Sorangiumaov<-aov(SorangiumplotPXG$pheno~SorangiumplotPXG$SCRI_RS_154747) 
summary(Sorangiumaov)
TukeyHSD(Sorangiumaov)

##Rhodanobacter
RhodanobacterplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=4, jitter=1, infer=T)
Rhodanobacteraov<-aov(RhodanobacterplotPXG$pheno~RhodanobacterplotPXG$SCRI_RS_154747) 
summary(Rhodanobacteraov) 
TukeyHSD(Rhodanobacteraov)
##Holophaga
HolophagaplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=12, jitter=1, infer=T)
Holophagaaov<-aov(HolophagaplotPXG$pheno~HolophagaplotPXG$SCRI_RS_154747) 
summary(Holophagaaov)
TukeyHSD(Holophagaaov)
##Tahibacter
TahibacterplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=15, jitter=1, infer=T)
Tahibacteraov<-aov(TahibacterplotPXG$pheno~TahibacterplotPXG$SCRI_RS_154747) 
summary(Tahibacteraov) 
TukeyHSD(Tahibacteraov)
##########
##Microbacterium
MicrobacteriumplotPXG<-plotPXG(wild, "SCRI_RS_189371", pheno.col=10, jitter=1, infer=T)
Microbacteriumaov<-aov(MicrobacteriumplotPXG$pheno~MicrobacteriumplotPXG$SCRI_RS_189371)
summary(Microbacteriumaov)
TukeyHSD(Microbacteriumaov)
##Streptomyces
StreptomycesplotPXG<-plotPXG(wild, "SCRI_RS_120720", pheno.col=18, jitter=1, infer=T)
Streptomycesaov<-aov(StreptomycesplotPXG$pheno~StreptomycesplotPXG$SCRI_RS_120720) 
summary(Streptomycesaov)
TukeyHSD(Streptomycesaov)

############################################################################################################################
##End

#####################################################################################################################
##QTL mapping Family level
#####################################################################################################################
###Read the file and define the population type
library(qtl)
library(devtools)

##Read.cross objects:https://www.rdocumentation.org/packages/qtl/versions/1.46-2/topics/read.cross

wild <- read.cross("csv", ".", "JH07_merged_HID144-Barke_Silva138_Family_enriched_Mapping.csv", BC.gen = 1, F.gen = 3)
###As the map position comes from genetic map wherein several markers map to same location, we might get error message but that's fine.
wild <- jittermap(wild)
summary(wild)
plot.map(wild)

###calculate the genotypic probability with almost default parameters
wild1<-calc.genoprob(wild, step=2.5,error.prob=0.01)
###similarity of the genotype
wild2<-sim.geno(wild, step=2.5, n.draws=64, error.prob=0.01)
###Different methods of analysis see the manual for details, specify the phenotype column
#Change here the number of phenotypes
out.em1 <- scanone(wild1, pheno.col =1:17) #number of phenotypes or taxa = 17
###get the summary of the analysis
summary(out.em1)

###Another method
out.hk1 <- scanone(wild1, pheno.col =1:17, method="hk")
###get the summary of the analysis
summary(out.hk1)
###Another method
out.imp1 <- scanone(wild2, pheno.col = 1:17, method="imp")
###get the summary of the analysis
summary(out.imp1)

#write.csv(out.em1, "JH20_Family_segregant_analysis_9K_out_em1_0721.csv")
#write.csv(out.hk1, "JH20_Family_segregant_analysis_9K_out_hk1_0721.csv")
#write.csv(out.imp1, "JH20_Family_segregant_analysis_9K_out_imp1_0721.csv")
#output LOD thresholds 
###################################################################################
###to get a LOD genome-wide significance threshold or genome-scan-adjusted p-values
####################################################################################
perm <- scanone(wild1, method="em", n.perm=1000,pheno.col = 1:17)
plot(perm)


#Significance thresholds may be obtained via the summary() function:
summary(perm)
summary(perm, alpha=c( 0.2))

#Most importantly, the permutation results may be used along with the scanone() results to have significance thresholds
#and p-values calculated automatically:
summary(out.em1, perms=perm, alpha=0.2, pvalues=TRUE)

#Significance thresholds were calculated by permutation, with thresholds of 20% and 5% genome wide significance 

#################################
##3. Test for epistasis
##################################

out22.em <-scantwo(wild1)

plot(out22.em)

summary(out22.em)


operm2 <-scantwo(wild1,method="hk",n.perm=100)

summary(operm2, alpha=c(0.05,0.20))


summary(out22.em)#NO evidence for epistasis

plot(out2.em, chr=c(3,5))

max(out22.em)


########################################################
##Confidence intervals
#######################################################
#To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL on chromosome of interest, type:

###Comamonaceae
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=6, expandtomarkers=FALSE)

###Polyangiaceae
bayesint(out.em1, chr="1H", out.2qtl, prob=0.95, lodcolumn=12, expandtomarkers=FALSE)
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=12, expandtomarkers=FALSE)

###Holophagaceae
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=10, expandtomarkers=FALSE)

###Rhodanobacteraceae
bayesint(out.em1, chr="3H", out.2qtl, prob=0.95, lodcolumn=13, expandtomarkers=FALSE)

##Sphingobacteriaceae
bayesint(out.em1, chr="5H", out.2qtl, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

###Streptomycetaceae
bayesint(out.em1, chr="7H", out.2qtl, prob=0.95, lodcolumn=15, expandtomarkers=FALSE)

###################################
##4. Fitting the model
###################################


##By taxa
######################
##Comamonaceae
qtl_Comamonaceae <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.Comamonaceae <- fitqtl(wild2, pheno.col=6, qtl=qtl_Comamonaceae , formula=y~Q1+Q2,method="imp", get.ests=TRUE))

##Polyangiaceae
qtl_Polyangiaceae <- makeqtl(wild2, chr=c("1H","3H"), pos=c(97.9,38.75))
summary(out.Polyangiaceae <- fitqtl(wild2, pheno.col=12, qtl=qtl_Polyangiaceae , formula=y~Q1+Q2,method="hk", get.ests=TRUE))

##Holophagaceae
qtl_Holophagaceae <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.Holophagaceae <- fitqtl(wild2, pheno.col=10, qtl=qtl_Holophagaceae , formula=y~Q1,method="hk", get.ests=TRUE))

##Rhodanobacteraceae
qtl_Rhodanobacteraceae <- makeqtl(wild2, chr=c("3H"), pos=c(38.75))
summary(out.Rhodanobacteraceae <- fitqtl(wild2, pheno.col=13, qtl=qtl_Rhodanobacteraceae , formula=y~Q1,method="imp", get.ests=TRUE))

####Sphingobacteriaceae
qtl_Sphingobacteriaceae <- makeqtl(wild2, chr=c("5H"), pos=c(158.3))
summary(out.Rhodanobacteraceae <- fitqtl(wild2, pheno.col=1, qtl=qtl_Sphingobacteriaceae , formula=y~Q1,method="imp", get.ests=TRUE))

##Streptomycetaceae
qtl_Streptomycetaceae <- makeqtl(wild2, chr=c("7H"), pos=c(133.9))
summary(out.2qtl <- fitqtl(wild2, pheno.col=15, qtl=qtl_Streptomycetaceae , formula=y~Q1,method="hk", get.ests=TRUE))


###########################################################################################################
##Marker regression: Allele composition at certain marker location:Plot the phenotype values versus the genotypes at a marker 
###########################################################################################################
#Comamonaceae
ComamonaceaeplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=6, jitter=1, infer=T)#,
Comamonaceaeaov<-aov(ComamonaceaeplotPXG$pheno~ComamonaceaeplotPXG$SCRI_RS_154747) 
summary(Comamonaceaeaov)  


##Polyangiaceae
PolyangiaceaeplotPXG<-plotPXG(wild, "BOPA1_1497_628", pheno.col=12, jitter=1, infer=T)#,
Polyangiaceaeaov<-aov(PolyangiaceaeplotPXG$pheno~PolyangiaceaeplotPXG$BOPA1_1497_628) 
summary(Polyangiaceaeaov)  

PolyangiaceaeplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=12, jitter=1, infer=T)#,
Polyangiaceaeaov<-aov(PolyangiaceaeplotPXG$pheno~PolyangiaceaeplotPXG$SCRI_RS_154747) 
summary(Polyangiaceaeaov)  

##Holophagaceae
HolophagaceaeplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=10, jitter=1, infer=T)
Holophagaceaeaov<-aov(HolophagaceaeplotPXG$pheno~HolophagaceaeplotPXG$SCRI_RS_154747) 
summary(Holophagaceaeaov)


##Rhodanobacteraceae
RhodanobacteraceaeplotPXG<-plotPXG(wild, "SCRI_RS_154747", pheno.col=13, jitter=1, infer=T)#,
Rhodanobacteraceaeaov<-aov(RhodanobacteraceaeplotPXG$pheno~RhodanobacteraceaeplotPXG$SCRI_RS_154747) 
summary(Rhodanobacteraceaeaov)  
TukeyHSD(Rhodanobacteraceaeaov) 
##Sphingobacteriaceae
SphingobacteriaceaeplotPXG<-plotPXG(wild, "BOPA1_6736_452", pheno.col=1, jitter=1, infer=T)#,
Sphingobacteriaceaeaov<-aov(SphingobacteriaceaeplotPXG$pheno~SphingobacteriaceaeplotPXG$BOPA1_6736_452) 
summary(Sphingobacteriaceaeaov) 
TukeyHSD(Sphingobacteriaceaeaov)
#Streptomycetaceae
StreptomycetaceaeplotPXG<-plotPXG(wild, "SCRI_RS_120015", pheno.col=15, jitter=1, infer=T)
Streptomycetaceaeaov<-aov(StreptomycetaceaeplotPXG$pheno~StreptomycetaceaeplotPXG$SCRI_RS_120015) 
summary(Streptomycetaceaeaov)
TukeyHSD(Streptomycetaceaeaov)
############################################################################################################################
##End