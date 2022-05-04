#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in: https://www.biorxiv.org/content/10.1101/2021.12.20.472907v1
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
###########################################################
library (pastecs)
library (ggplot2)
library(dplyr)
library(MASS)
library(ggfortify)
library(ggthemes)
library(colorspace)
library(grDevices)
library(devtools)
##########################
#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

setwd ("")
getwd()

#################################################
################################################
##Supplementary Table 4 Root_architecture
###############################################
##############################################

#Import data 

Root_arq_CE3 <- read.delim("ROOT_ARCH_CE_3_Entire_root_F15.txt", sep = "\t", stringsAsFactors = TRUE) 

by(Root_arq_CE3, Root_arq_CE3$Genotype, summary)

by(Root_arq_CE3, Root_arq_CE3$Genotype, stat.desc)


########################## Shoot_weight(g)
#############################################

shapiro.test(Root_arq_CE3$Shoot.weight)

qqnorm (Root_arq_CE3$Shoot.weight) + qqline(Root_arq_CE3$Shoot.weight)
hist(Root_arq_CE3$Shoot.weight)

Shoot.weight_aov<-aov(Shoot.weight~Genotype, data=Root_arq_CE3)
summary(Shoot.weight_aov)

########################## Root.weight(g)
############################################

shapiro.test(Root_arq_CE3$Root.weight)
hist(Root_arq_CE3$Root.weight)

kruskal.test(Root.weight~Genotype , data=Root_arq_CE3)

########################## Root Length (cm)
##########################################

shapiro.test(Root_arq_CE3$RootLength)
hist(Root_arq_CE3$RootLength)

Root_length_aov<-aov(RootLength~Genotype, data=Root_arq_CE3)
summary(Root_length_aov)

##########################  Primary Root Length PRL (cm)
#######################################################
shapiro.test(Root_arq_CE3$PRL)
hist(Root_arq_CE3$PRL)

kruskal.test(PRL~Genotype , data=Root_arq_CE3)

######################## Surface Area (cm2)
##########################################

shapiro.test(Root_arq_CE3$SurfArea)
hist(Root_arq_CE3$SurfArea)

SurfArea_aov<-aov(SurfArea~Genotype, data=Root_arq_CE3)
summary(SurfArea_aov)

########################### Average Diameter (cm)
###########################################

shapiro.test(Root_arq_CE3$AvgDiam)
hist(Root_arq_CE3$AvgDiam )

AvgDiam_aov<-aov(AvgDiam~Genotype, data=Root_arq_CE3)
summary(AvgDiam_aov)

######################## Root Volume (cm3)
##############################################

shapiro.test(Root_arq_CE3$RootVol)
hist(Root_arq_CE3$RootVol)

RootVol_aov<-aov(RootVol~Genotype, data=Root_arq_CE3)
summary(RootVol_aov)

############################# Tips
######################################

shapiro.test(Root_arq_CE3$Tips)
hist(Root_arq_CE3$Tips)

Tips_aov<-aov(Tips~Genotype, data=Root_arq_CE3)
summary(Tips_aov)

############################ Forks
###################################

shapiro.test(Root_arq_CE3$Forks)
hist(Root_arq_CE3$Forks)

Forks_aov<-aov(Forks~Genotype, data=Root_arq_CE3)
summary(Forks_aov)

################################ Crossings
###########################################

shapiro.test(Root_arq_CE3$Crossings)
hist(Root_arq_CE3$Crossings)

Crossings_aov<-aov(Crossings~Genotype, data=Root_arq_CE3)
summary(Crossings_aov)

############################ Specific Root length
###################################################

shapiro.test(Root_arq_CE3$SRL)
hist(Root_arq_CE3$SRL)

SRL_aov<-aov(SRL~Genotype, data=Root_arq_CE3)
summary(SRL_aov)

############################Root Density
###########################################

shapiro.test(Root_arq_CE3$RootDen)
hist(Root_arq_CE3$RootDen)    

kruskal.test(RootDen~Genotype , data=Root_arq_CE3)

##End
############################################################################################################################################