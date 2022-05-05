#Ref to the ARTICLE 
# 
#  Code to compute calculations presented in:https://www.biorxiv.org/content/10.1101/2021.12.20.472907v3
#  Supplementary Figure 10, Database 2
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
library (pastecs) 
library(dplyr)
library(FSA)
##########################
#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

#import data 

Root_exu <- read.delim("Polar_Ery_F15_0322.txt", sep = "\t", stringsAsFactors = TRUE)

by(Root_exu, Root_exu$Genotype, summary)
by(Root_exu, Root_exu$Genotype, stat.desc)#library pastecs

##############################################################################
##Polar compounds
##############################################################################
##############################################################################
#L-valine
shapiro.test(Root_exu$L.valine.TMS2.)

test<-aov(Root_exu$L.valine.TMS2.~ Root_exu$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test

#L-Leucine
shapiro.test(Root_exu$L.leucine.TMS2.)

kruskal.test(Root_exu$L.leucine.TMS2.~ Root_exu$Genotype)

dunnTest(L.leucine.TMS2.~ Genotype,
         data=Root_exu,method="bh")

##Glycerol
shapiro.test(Root_exu$Glycerol.TMS3.)

kruskal.test(Root_exu$Glycerol.TMS3.~ Root_exu$Genotype)

dunnTest(P_1263_Glycerol.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

#L-isoleucine
shapiro.test(Root_exu$L.isoleucine.TMS2.)

kruskal.test(Root_exu$L.isoleucine.TMS2.~ Root_exu$Genotype)

dunnTest(L.isoleucine.TMS2. ~ Genotype,
         data=Root_exu,method="bh")


#L-proline
shapiro.test(Root_exu$L.proline.TMS2.)

test<-aov(Root_exu$L.proline.TMS2.~ Root_exu$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test

#Succinic acid
shapiro.test(Root_exu$Succinic.acid)

kruskal.test(Root_exu$Succinic.acid~ Root_exu$Genotype)

dunnTest(Succinic.acid ~ Genotype,
         data=Root_exu,method="bh")

#L-serine
shapiro.test(Root_exu$L.serine.TMS3.)

kruskal.test(Root_exu$L.serine.TMS3.~ Root_exu$Genotype)

dunnTest(L.serine.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

#L-threonine
shapiro.test(Root_exu$L.threonine.TMS3.)

test<-aov(Root_exu$L.threonine.TMS3.~ Root_exu$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test

#Malic_acid

shapiro.test(Root_exu$Malic_acid.TMS3.)

kruskal.test(Root_exu$Malic_acid.TMS3.~ Root_exu$Genotype)

dunnTest(Malic_acid.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

#L-proline
shapiro.test(Root_exu$X5.oxo.proline.TMS2.)

kruskal.test(Root_exu$X5.oxo.proline.TMS2.~ Root_exu$Genotype)

dunnTest(X5.oxo.proline.TMS2. ~ Genotype,
         data=Root_exu,method="bh")

#L-aspartic acid

shapiro.test(Root_exu$L.aspartic_acid.TMS3.)
kruskal.test(Root_exu$L.aspartic_acid.TMS3.~ Root_exu$Genotype)

dunnTest(L.aspartic_acid.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

# Gamma-aminobutyric acid
shapiro.test(Root_exu$gamma.aminobutyric_acid.TMS3.)

kruskal.test(Root_exu$gamma.aminobutyric_acid.TMS3.~ Root_exu$Genotype)

dunnTest(gamma.aminobutyric_acid.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

#L.glutamic_acid 
shapiro.test(Root_exu$ L.glutamic_acid)

kruskal.test(Root_exu$L.glutamic_acid~ Root_exu$Genotype)

dunnTest(L.glutamic_acid ~ Genotype,
         data=Root_exu,method="bh")

#Hydroxyphenylacetic.acid
shapiro.test(Root_exu$X4.hydroxyphenylacetic.acid.TMS2.)

kruskal.test(Root_exu$X4.hydroxyphenylacetic.acid.TMS2.~ Root_exu$Genotype)

dunnTest(X4.hydroxyphenylacetic.acid.TMS2. ~ Genotype,
         data=Root_exu,method="bh")
#L.asparagine

shapiro.test(Root_exu$L.asparagine.TMS3.)

kruskal.test(Root_exu$L.asparagine.TMS3.~ Root_exu$Genotype)


#Ribose
shapiro.test(Root_exu$Ribose)

kruskal.test(Root_exu$Ribose~ Root_exu$Genotype)

dunnTest(Ribose ~ Genotype,
         data=Root_exu,method="bh")

#L.asparagine
shapiro.test(Root_exu$L.asparagine.TMS3.)

kruskal.test(Root_exu$L.asparagine.TMS3.~ Root_exu$Genotype)

dunnTest(L.asparagine.TMS3. ~ Genotype,
         data=Root_exu,method="bh")

#phenylpropanolamine
shapiro.test(Root_exu$phenylpropanolamine.TMS2.)

kruskal.test(Root_exu$phenylpropanolamine.TMS2.~ Root_exu$Genotype)

dunnTest(phenylpropanolamine.TMS2. ~ Genotype,
         data=Root_exu,method="bh")

#Putrescine

shapiro.test(Root_exu$Putrescine.TMS4.)

test<-aov(Root_exu$Putrescine.TMS4.~ Root_exu$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test

#L.glutamine
shapiro.test(Root_exu$L.glutamine.TMS3.)

kruskal.test(Root_exu$L.glutamine.TMS3.~ Root_exu$Genotype)

dunnTest(L.glutamine.TMS3. ~ Genotype,
         data=Root_exu,method="bh")
#Azelaic_acid

shapiro.test(Root_exu$Azelaic_acid.TMS2.)

kruskal.test(Root_exu$Azelaic_acid.TMS2.~ Root_exu$Genotype)

dunnTest(Azelaic_acid.TMS2. ~ Genotype,
         data=Root_exu,method="bh")

#Citric acid

shapiro.test(Root_exu$ Citric_acid.TMS4.)

kruskal.test(Root_exu$ Citric_acid.TMS4.~ Root_exu$Genotype)

dunnTest( Citric_acid.TMS4. ~ Genotype,
          data=Root_exu,method="bh")

#Unknown3

shapiro.test(Root_exu$Unknown3)

kruskal.test(Root_exu$Unknown3~ Root_exu$Genotype)

dunnTest(Unknown3 ~ Genotype,
         data=Root_exu,method="bh")

#Sugar dimer 1

shapiro.test(Root_exu$Sugar.dimer)

kruskal.test(Root_exu$Sugar.dimer~ Root_exu$Genotype)

dunnTest(Sugar.dimer~ Genotype,
         data=Root_exu,method="bh")


#Gluconic_acid
shapiro.test(Root_exu$Gluconic_acid.TMS6.)

kruskal.test(Root_exu$Gluconic_acid.TMS6.~ Root_exu$Genotype)

dunnTest(Gluconic_acid.TMS6.~ Genotype,
         data=Root_exu,method="bh")


#Gulose

shapiro.test(Root_exu$Gulose.TMS5.)

kruskal.test(Root_exu$Gulose.TMS5.~ Root_exu$Genotype)

dunnTest(Gulose.TMS5.~ Genotype,
         data=Root_exu,method="bh")

#Inositol

shapiro.test(Root_exu$Inositol.TMS6.)

kruskal.test(Root_exu$Inositol.TMS6.~ Root_exu$Genotype)

dunnTest(Inositol.TMS6.~ Genotype,
         data=Root_exu,method="bh")

#glycerol.galactoside
shapiro.test(Root_exu$X2.O.glycerol.galactoside.TMS6.)

kruskal.test(Root_exu$X2.O.glycerol.galactoside.TMS6.~ Root_exu$Genotype)

dunnTest(X2.O.glycerol.galactoside.TMS6.~ Genotype,
         data=Root_exu,method="bh")

#Glucuronic_acid
shapiro.test(Root_exu$Glucuronic_acid.TMS5.)

test<-aov(Root_exu$Glucuronic_acid.TMS5.~ Root_exu$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test

#Sucrose
shapiro.test(Root_exu$Sucrose.TMS8.)

kruskal.test(Root_exu$Sucrose.TMS8.~ Root_exu$Genotype)

dunnTest(Sucrose.TMS8.~ Genotype,
         data=Root_exu,method="bh")
#Unknown4 
shapiro.test(Root_exu$Unknown4 )

kruskal.test(Root_exu$Unknown4 ~ Root_exu$Genotype)

dunnTest(Unknown4 ~ Genotype,
         data=Root_exu,method="bh")

#Unknown5
shapiro.test(Root_exu$Unknown5 )

kruskal.test(Root_exu$Unknown5 ~ Root_exu$Genotype)

dunnTest(Unknown5 ~ Genotype,
         data=Root_exu,method="bh")

#Maltose
shapiro.test(Root_exu$Maltose_MEOX.TMS8.)

kruskal.test(Root_exu$Maltose_MEOX.TMS8.~ Root_exu$Genotype)

dunnTest(Maltose_MEOX.TMS8.~ Genotype,
         data=Root_exu,method="bh")

#Melibiose

shapiro.test(Root_exu$Melibiose.TMS8.)

kruskal.test(Root_exu$Melibiose.TMS8.~ Root_exu$Genotype)

dunnTest(Melibiose.TMS8.~ Genotype,
         data=Root_exu,method="bh")

#Raffinose

shapiro.test(Root_exu$Raffinose.TMS8.)

kruskal.test(Root_exu$Raffinose.TMS8.~ Root_exu$Genotype)

dunnTest(Raffinose.TMS8.~ Genotype,
         data=Root_exu,method="bh")

#Fructose

shapiro.test(Root_exu$Frutose_total)

kruskal.test(Root_exu$Frutose_total~ Root_exu$Genotype)

dunnTest(Frutose_total~ Genotype,
         data=Root_exu,method="bh")

#Glucose

shapiro.test(Root_exu$Glucose_total)

kruskal.test(Root_exu$Glucose_total~ Root_exu$Genotype)

dunnTest(Glucose_total~ Genotype,
         data=Root_exu,method="bh")

##############################################################################
##Non- polar compounds
##############################################################################
##############################################################################
Root_exu <- read.delim("NonPolar_0322_2.txt", sep = "\t", stringsAsFactors = TRUE)

by(Root_exu, Root_exu$Genotype, summary)
by(Root_exu, Root_exu$Genotype, stat.desc)#library pastecs

##NP_1243_Unknown2 
shapiro.test(F15$NP_1243_Unknown2)

kruskal.test(F15$NP_1243_Unknown2~ F15$Genotype)

#NP_1485_Unknown4

shapiro.test(F15$NP_1485_Unknown4)

kruskal.test(F15$NP_1485_Unknown4~ F15$Genotype)

dunnTest(NP_1485_Unknown4~ Genotype,
         data=F15,method="bh")

#Methyl.tetradecanoic_acid
shapiro.test(F15$NP_1807_12.methyl.tetradecanoic_acid)

test<-aov(F15$NP_1807_12.methyl.tetradecanoic_acid~ F15$Genotype) # print results of anova
summary(test)

tukey.test <- TukeyHSD(test)
tukey.test


#NP_1865_Unknown5

shapiro.test(F15$NP_1865_Unknown5)

kruskal.test(F15$NP_1865_Unknown5~ F15$Genotype)

dunnTest(NP_1865_Unknown5~ Genotype,
         data=F15,method="bh")

#NP_1897_methyl.pentadecanoic_acid_ME

shapiro.test(F15$NP_1897_methyl.pentadecanoic_acid_ME)

test<-aov(F15$NP_1897_methyl.pentadecanoic_acid_ME~ F15$Genotype) # print results of anova
summary(test)


#Hexadecanoic_acid

shapiro.test(F15$Hexadecenoic.acid)

kruskal.test(F15$Hexadecenoic.acid~ F15$Genotype)


##NP_2008_cyclopropaneoctanoic_acid_2.hex
shapiro.test(F15$NP_2008_cyclopropaneoctanoic_acid_2.hex)

kruskal.test(F15$NP_2008_cyclopropaneoctanoic_acid_2.hex~ F15$Genotype)

dunnTest(NP_2008_cyclopropaneoctanoic_acid_2.hex~ Genotype,
         data=F15,method="bh")


####9,12-octadecadienoic acid

shapiro.test(F15$NP_2097_9.12.octadecadienoic_acid_ME)

test<-aov(F15$NP_2097_9.12.octadecadienoic_acid_ME~ F15$Genotype) # print results of anova
summary(test)


##Octadecatrienoic_acid

shapiro.test(F15$Octadecatrienoic_acid)

test<-aov(F15$Octadecatrienoic_acid~ F15$Genotype) # print results of anova
summary(test)


#NP_2186_3.7.11.15.tetramethyl.hexadecanol.

shapiro.test(F15$NP_2186_3.7.11.15.tetramethyl.hexadecanol)

kruskal.test(F15$NP_2186_3.7.11.15.tetramethyl.hexadecanol~ F15$Genotype)



#NP_2336_9.12.15.Octadecatrienoic_acid.gly

shapiro.test(F15$NP_2336_9.12.15.Octadecatrienoic_acid.gly)

kruskal.test(F15$NP_2336_9.12.15.Octadecatrienoic_acid.gly~ F15$Genotype)


#NP_2893_Unknown10

shapiro.test(F15$NP_2893_Unknown10)

kruskal.test(F15$NP_2893_Unknown10~ F15$Genotype)

dunnTest(NP_2893_Unknown10~ Genotype,
         data=F15,method="bh")

#Methyl.cholesterol

shapiro.test(F15$NP_3149_Methyl.cholesterol)

kruskal.test(F15$NP_3149_Methyl.cholesterol~ F15$Genotype)


#Stigmasterol

shapiro.test(F15$NP_3187_Stigmasterol)

kruskal.test(F15$NP_3187_Stigmasterol~ F15$Genotype)


#Alpha.sitosterol

shapiro.test(F15$NP_3269_alpha.sitosterol)

kruskal.test(F15$NP_3269_alpha.sitosterol~ F15$Genotype)

########################################################################################
##End
