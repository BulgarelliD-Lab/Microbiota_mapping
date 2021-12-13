#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:
#  Supplementary Fig.7
#  Revision 11/21 
#  c.m.z.escuderomartinez@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

library ("ggplot2")
library("dplyr")
library("MASS")
library("ggfortify")
library("ggthemes")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMRplus")

#retrieve R and package versions for the reproducibility of the code
sessionInfo()

#################################################
################################################
##Supplementary Fig. 7 Elemental analysis
###############################################
##############################################

#Import data file all lines
EA_results <- read.delim("Elemental_analyzer_F15_results.txt", sep = "\t")

##################################
##Subset by timepoint
##################################

AE_3weeks<-subset(EA_results,TP == "3")
AE_2weeks<-subset(EA_results,TP=="2")

##################################################################
###Per genotype and timepoint
####################################################################

## Nitrogen 3 weeks 
############################################################
shapiro.test(AE_3weeks$N)


kruskal.test(N~Genotype, data=AE_3weeks)


KDn3<-posthoc.kruskal.dunn.test (x= AE_3weeks$N, g= AE_3weeks$Genotype, p.adjust.method="none")

## Carbon 3 weeks 
############################################################
shapiro.test(AE_3weeks$C)

kruskal.test(C~Genotype, data=AE_3weeks)


KDc3<-posthoc.kruskal.dunn.test (x= AE_3weeks$C, g= AE_3weeks$Genotype, p.adjust.method="none")

KDc3

##Plotting
###############################################

#Order the levels according to a defined order 
AE_3weeks$Genotype <- ordered(AE_3weeks$Genotype, levels=c("NC", "Barke", "Int 17","Int 52"))  

p<-ggplot(AE_3weeks, aes(x=Genotype, y=N, fill=Genotype))+ 
  geom_bar(position=position_dodge(0.8),stat='identity')+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,20)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

p+scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)


p<-ggplot(AE_3weeks, aes(x=Genotype, y=C, fill=Genotype))+ 
  geom_bar(position=position_dodge(0.8),stat='identity')+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,30)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

p+scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)
###################################################################################


## Nitrogen 2 weeks 
############################################################
shapiro.test(AE_2weeks$N)

N_2weeks_aov<-aov(N~Genotype, data=AE_2weeks)
summary(N_2weeks_aov)

## Carbon 2 weeks 
############################################################
shapiro.test(AE_2weeks$C)

kruskal.test(C~Genotype, data=AE_2weeks)


KDc2<-posthoc.kruskal.dunn.test (x= AE_2weeks$C, g= AE_2weeks$Genotype, p.adjust.method="none")

KDc2

##Plotting
###############################################

##Order the levels according to a defined order 
AE_2weeks$Genotype <- ordered(AE_2weeks$Genotype, levels=c("NC", "Barke", "Int 17","Int 52"))#,"Int 19","Int 56"))  

p<-ggplot(AE_2weeks, aes(x=Genotype, y=N, fill=Genotype))+ 
  geom_bar(position=position_dodge(0.8),stat='identity')+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,20)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


p+scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

p<-ggplot(AE_2weeks, aes(x=Genotype, y=C, fill=Genotype))+ 
  geom_bar(position=position_dodge(0.8),stat='identity')+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,30)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

p+scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

##End
#############################################################################################################################################








