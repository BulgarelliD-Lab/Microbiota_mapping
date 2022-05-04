#####################################################################################
#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:https://www.biorxiv.org/content/10.1101/2021.12.20.472907v1
#  SFigure 9
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

library ("ggplot2")
library("dplyr")
library("MASS")
library("ggfortify")
library("ggthemes")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMR")
library("lsmeans")
library("dunn.test")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

##working directory

setwd ("")
getwd()

#####################################################
####################################################
##Supplementary Fig 9 yield 
###################################################
##################################################

#Import yield data file all lines
HIF_pheno <- read.delim("HIF_all4_reps_Int_1121.txt", sep = "\t") 

HIF_pheno<-subset(HIF_pheno, Genotype == "Barke"| Genotype == "Int_17"|Genotype =="Int_52")

###############################
##TGW
####################################
#plotting TGW

HIF_pheno$Genotype <- ordered(HIF_pheno$Genotype, levels=c("Barke", "Int_17","Int_52"))

p<-ggplot(HIF_pheno, aes(x=Genotype, y=TGW, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(25,60)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
p+facet_grid(~Rep)+scale_fill_manual(values=c("#0072B2", "#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

#Stast 
shapiro.test(HIF_pheno$TGW)

TGW_aov<-aov(TGW~Genotype*Rep, data=HIF_pheno)
summary(TGW_aov)

TukeyHSD(TGW_aov)

R1<-subset(HIF_pheno, Rep == "1R")
TGW_R1_aov<-aov(TGW~Genotype, data=R1)
summary(TGW_R1_aov)

R2<-subset(HIF_pheno, Rep == "2R")
TGW_R2_aov<-aov(TGW~Genotype, data=R2)
summary(TGW_R2_aov)
TukeyHSD(TGW_R2_aov)

R3<-subset(HIF_pheno, Rep == "3R")
TGW_R3_aov<-aov(TGW~Genotype, data=R3)
summary(TGW_R3_aov)
TukeyHSD(TGW_R3_aov)

R4<-subset(HIF_pheno, Rep == "4R")
TGW_R4_aov<-aov(TGW~Genotype, data=R4)
summary(TGW_R4_aov)
TukeyHSD(TGW_R4_aov)


###############################
##Weight
####################################
#plotting total weight

##Order the levels according to a defined order 
HIF_pheno$Genotype <- ordered(HIF_pheno$Genotype, levels=c("Barke", "Int_17","Int_52"))

p<-ggplot(HIF_pheno, aes(x=Genotype, y=Weight, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,1.5)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

p+facet_grid(~Rep)+scale_fill_manual(values=c("#0072B2", "#56B4E9","#E69F00","#D55E00","#D55E00"))+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)


  geom_point(aes(fill=Genotype, shape = Rep, group = Genotype),size=5,position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)+
  scale_shape_manual(values=c(21,22,23))

#stast 
shapiro.test(HIF_pheno$Weight)


Weight_aov<-aov(Weight~Genotype*Rep, data=HIF_pheno)
summary(Weight_aov)

TukeyHSD(Weight_aov)

R1<-subset(HIF_pheno, Rep == "1R")
Weight_R1_aov<-aov(Weight~Genotype, data=R1)
summary(Weight_R1_aov)

R2<-subset(HIF_pheno, Rep == "2R")
Weight_R2_aov<-aov(Weight~Genotype, data=R2)
summary(Weight_R2_aov)
TukeyHSD(Weight_R2_aov)

R3<-subset(HIF_pheno, Rep == "3R")
Weight_R3_aov<-aov(Weight~Genotype, data=R3)
summary(Weight_R3_aov)
TukeyHSD(Weight_R3_aov)

R4<-subset(HIF_pheno, Rep == "4R")
Weight_R4_aov<-aov(Weight~Genotype, data=R4)
summary(Weight_R4_aov)
TukeyHSD(Weight_R4_aov)

##End
##############################################################################################################################################
