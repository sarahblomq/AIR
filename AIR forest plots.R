
library(tidyverse)
library(readxl)
library(effects)
library(RNOmni)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr)
library(kinship2)
library(coxme)
library(MuMIn)

setwd("~/Chilton Laboratory/Gene-Diet Cohorts/AIR/Final AIR Datasheets")
dgf<-read_excel("AIR FA geno conc final.xlsx", na = "NA")
dp <- read_excel("AIR phenotypes - vars.xlsx", na = "NA")
dped <- read_excel ("AIR_Pedigree_forAnalysis.xlsx")
DP <- merge(dp, dped, by = "id")
d <- merge.data.frame(dgf, DP)

#set rs ids & sex as factor
#genotype as integer
d$sex <-  as.factor(d$sex)


#currently coded by ANCESTRAL allele, except rs174538 which is flipped because of 1000G haplotypes

d$rs1535 <-  factor(d$rs1535, levels = c("A", "AG", "G"))
d$rs174455 <-  factor(d$rs174455, levels = c("A", "AG", "G"))
d$rs174537 <-  factor(d$rs174537, levels = c("G", "GT", "T"))
d$rs174538 <-  factor(d$rs174538, levels = c("G", "GA", "A"), labels = c("G", "GA", "A"))
d$rs174546 <-  factor(d$rs174546, levels = c("C", "CT", "T"))
d$rs174547 <-  factor(d$rs174547, levels = c("T", "CT", "C"), labels = c("T", "TC", "C"))
d$rs174548 <-  factor(d$rs174548, levels = c("C", "GC", "G"), labels = c("C", "CG", "G"))
d$rs174554 <-  factor(d$rs174554, levels = c("A", "AG", "G"))
d$rs174570 <-  factor(d$rs174570, levels = c("T", "CT", "C"), labels = c("T", "TC", "C"))
d$rs174576 <-  factor(d$rs174576, levels = c("C", "CA", "A"), labels = c("C", "CA", "A"))
d$rs174579 <-  factor(d$rs174579, levels = c("T", "TC", "C"), labels = c("T", "TC", "C"))
d$rs174594 <-  factor(d$rs174594, levels = c("A", "CA", "C"), labels = c("A", "AC", "C"))
d$rs174602 <-  factor(d$rs174602, levels = c("T", "TC", "C"), labels = c("T", "CT", "C"))
d$rs7115739 <-  factor(d$rs7115739, levels = c("T", "GT", "G"), labels = c("T", "TG", "G"))

#ELOVL2
#According to LDhap/1000G PEL main hap, rs3734398 = C, rs3798713 = C, rs953413 = A, rs3796719 = C
d$rs3734398 <-  factor(d$rs3734398, levels = c("T", "CT", "C"), labels = c("T", "TC", "C"))
d$rs3798713 <-  factor(d$rs3798713, levels = c("G", "CG", "C"), labels = c("G", "GC", "C"))
d$rs3798719 <-  factor(d$rs3798719, levels = c("T", "TC", "C"), labels = c("T", "CT", "C"))
d$rs953413 <-  factor(d$rs953413, levels = c("G", "AG", "A"), labels = c("G", "GA", "A"))
d$rs1570069 <-  factor(d$rs1570069, levels = c("A", "AG", "G"))
d$rs2281591 <-  factor(d$rs2281591, levels = c("G", "AG", "A"), labels = c("G", "GA", "A"))

#ELOVL5
#According to LDhap/1000G PEL main hap, rs7744440  = T, rs9357760 = A, rs2397142 = C
d$rs2397142 <-  factor(d$rs2397142, levels = c("G", "CG", "C"), labels = c("G", "GC", "C"))
d$rs9357760 <-  factor(d$rs9357760, levels = c ("G", "GA", "A"), labels = c("G", "GA", "A"))
d$rs7744440 <-  factor(d$rs7744440, levels = c("G", "GT", "T"), labels = c("G", "GT", "T"))

#genotype as integer for additive model

d$genors174455 <-  as.integer(d$rs174455)-1
d$genors7115739 <-  as.integer(d$rs7115739)-1
d$genors174602 <-  as.integer(d$rs174602)-1
d$genors174594 <-  as.integer(d$rs174594)-1
d$genors174579 <-  as.integer(d$rs174579)-1
d$genors174576 <-  as.integer(d$rs174576)-1
d$genors1535 <-  as.integer(d$rs1535)-1
d$genors174570 <-  as.integer(d$rs174570)-1
d$genors174554 <-  as.integer(d$rs174554)-1
d$genors174548 <-  as.integer(d$rs174548)-1
d$genors174547 <-  as.integer(d$rs174547)-1
d$genors174546 <-  as.integer(d$rs174546)-1
d$genors174538 <-  as.integer(d$rs174538)-1
d$genors174537 <-  as.integer(d$rs174537)-1
d$genors3798719 <-  as.integer(d$rs3798719)-1
d$genors1570069 <-  as.integer(d$rs1570069)-1
d$genors953413 <-  as.integer(d$rs953413)-1
d$genors3798713 <-  as.integer(d$rs3798713)-1
d$genors2281591 <-  as.integer(d$rs2281591)-1
d$genors3734398 <-  as.integer(d$rs3734398)-1
d$genors2397142 <-  as.integer(d$rs2397142)-1
d$genors9357760 <-  as.integer(d$rs9357760)-1
d$genors7744440 <-  as.integer(d$rs7744440)-1

#checking out the data & transforming if necessary
#ski did not want FA data log transformed as to keep the % total area representative

#log transform for normalization and unit consistency
phenoList = names(d)[c(4, 6:33, 35:38, 41:43, 81:102)]
phenoListLOG <- phenoList[c(1:18, 19:23, 24:29, 30:36, 37:58)]
phenoLL2 <- phenoList[c(15)]
d <- d %>% mutate(across(c(phenoListLOG), log))

#C18_2 needs an additional transformation
d <- d %>% mutate(across(c(phenoLL2), log))

#kinship matrix for adjusting for relatedness

ped <- pedigree(id=dped$ego_id, dadid=dped$fa_id, momid=dped$mo_id, famid=dped$fam_id, sex=dped$pedsex)
kin <- makekinship(dped$fam_id, dped$ego_id, dped$fa_id, dped$mo_id)

#add new t2d variable to original data set AFTER log transformation
d <- d %>%
  mutate(t2d = case_when(hba1c < 1.871802 & fast_glu < 4.836282 & Glu_2hr < 5.298317 ~ 'non-T2D',
                         hba1c >= 1.740466 & hba1c <= 1.856298 |  fast_glu >= 4.60517 & fast_glu <= 4.828314 |
                           Glu_2hr >= 4.941642 & Glu_2hr <= 5.293305 ~ 'pre-T2D',
                         hba1c >= 1.871802 | fast_glu >= 4.836282 | Glu_2hr >= 5.298317 ~ 'T2D',
  ))

#for not log transformed data
d <- d %>%
  mutate(t2d = case_when(hba1c < 6.5 & fast_glu < 126 & Glu_2hr < 200 ~ 'non-T2D',
                         hba1c >= 5.7 & hba1c <= 6.4 |  fast_glu >= 100 & fast_glu <= 125 |
                           Glu_2hr >= 140 & Glu_2hr <= 199 ~ 'pre-T2D',
                         hba1c >= 6.5 | fast_glu >= 126 | Glu_2hr >= 200 ~ 'T2D',
  ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#loop for forest plot master spreadsheet

snpList = names(d)[c(108, 110:111, 113:114, 116:130)]
phenoList = names(d)[c(4, 6:11, 13:31, 33, 35:38, 41:43, 81:98, 100:102)]
phenoEstList = c("C12", "C14", "C14_1", "C15", "C15_1", "C16", "C16_1", "C17_1", "C18", "C18_1cis", "C18_1trans",
                 "C18_1vac", "C18_2trans", "C18_2cis", "C18_3n3", "C18_3n6", "C18_4n3", "C20", "C20_1", "C20_2", "C20_3n3",
                 "C20_3n6", "C20_4n3", "C20_4n6", "C20_5n3", "C21", "C22_1", "C22_4n6", "C22_5n3", "C22_5n6", "C22_6", "C24_1", "C4", "C6",
                 "fatmass", "wc", "hip_circum", "chol", "trg", "hdl", "ldl", "vldl", "sbp", "dbp", "alt", "ast",
                 "fast_ins", "fast_glu", "Glu_2hr", "matsuda_index", "homair", "dio", "adiponectin", "height", "weight")
phenoEstList2 = ("C18_2cis")
phenoUnitList <- rep(c("ng/ml", "% (real)", "cm", "mg/dl", "mmHg", "IU/L", "uIu/ml", "mg/dl", "index units", "ug/ml", "cm", "kg"), 
                     times=c(34, 1, 2, 5, 2, 2, 1, 2, 3, 1, 1, 1))

#note: coefficients will need to be transformed if log transformed previously
#for calculating z value: fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))
#for calculating std error: sqrt(diag(vcov(tmp.lm)))
#for calculating p value: 2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))))
#for calculating CIs, from Yann: est +- (1.96 * std error)

#M/F split to check anthropometric
#add/remove sex var from tmp.formula as necessary
MFP <- d[which(d$sex=="M"),]
FFP <- d[which(d$sex=="F"),]

FP<-NULL; tmp.vec<-NULL
for (i in 1:length(snpList)) {
  for (j in 1:length(phenoList)) {
    tmp.formula <- paste0(phenoList[j], " ~ age + sex + bmi + (1|ego_id) + t2d + ", snpList[i])
    tmp.lm <- lmekin(as.formula(tmp.formula), data=FFP)
    tmp.p <- as.numeric(2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm))))))[6]
    if (phenoList[j] == phenoEstList[j]) {
      SE <- (exp((as.numeric(sqrt(diag(vcov(tmp.lm)))[6])))-1)*100
      est <- ((exp(as.numeric(tmp.lm$coefficients$fixed[6])))-1)*100
      conf.low <- est - (1.96 * SE)
      conf.high <- est + (1.96 * SE)
      tmp.vec <- c(snpList[i], phenoList[j], est, "%", conf.low, conf.high, tmp.p)
    } else if (phenoList[j] == phenoEstList2){
      SE <- (exp(exp((as.numeric(sqrt(diag(vcov(tmp.lm)))[6]))))-1)*100
      estLA <- (exp(exp(as.numeric(tmp.lm$coefficients$fixed[6])))-1)*100
      conf.low <- estLA - (1.96 * SE)
      conf.high <- estLA + (1.96 * SE)
      tmp.vec <- c(snpList[i], phenoList[j], estLA, "%", conf.low, conf.high, tmp.p) 
    } else if (phenoEstList[j] == "NA"){
      SE <- as.numeric(sqrt(diag(vcov(tmp.lm)))[6])
      conf.low <- as.numeric(tmp.lm$coefficients$fixed[6]) - (1.96 * SE)
      conf.high <- as.numeric(tmp.lm$coefficients$fixed[6]) + (1.96 * SE)
      tmp.vec <- c(snpList[i], phenoList[j], as.numeric(tmp.lm$coefficients$fixed[6]), phenoUnitList[j], conf.low, conf.high, tmp.p)
    }
    print(snpList[i])
    print(phenoList[j])
    FP <- rbind(FP, tmp.vec)
  }
}

FP <- as.data.frame(FP)
names(FP) = c("SNP","Phenotype","Est", "Est unit", "Conf.low", "Conf.high", "p_val")
rownames(FP)<-NULL
View(FP)

write.csv(FP,"~/Chilton Laboratory/Gene-Diet Cohorts/AIR/byT2D/FP - forest plot master conc.csv", row.names = TRUE)
################################################################################################
#forest plot figs for additive model

phenoList = names(d)[c(4, 6:11, 13:31, 33, 35:38, 41:43, 81:98, 100:102)]
phenoNames = c("Dodecanoic acid (C12)", "Myristic Acid (C14)", "Tetradecenoic Acid (C14:1)", "Pentadecylic Acid (C15)",
               "Pentadecenoic Acid (C15:1)", "Palmitic Acid (C16)", "Palmitoleic Acid (C16:1)",
               "Margaric Acid (C17:1)", "Stearic Acid (C18)", "Oleic Acid (C18:1 cis)", "Elaidic Acid (C18:1 trans)",
               "Vaccenic Acid (C18:1)", "Linolelaid Acid (C18:2 trans)", "Linoleic Acid (C18:2 cis)", "Alpha Linolenic Acid (C18:3 n-3)",
               "y-Linolenic Acid (C18:3 n-6)", "Stearidonic acid (C18:4 n-3)", "Arachidic Acid (C20)", "Paullinic Acid (C20:1)",
               "Dihomo-Linolenic Acid (C20:2)", "Dihomo-a-Linolenic Acid (C20:3 n-3)", "Dihomo-γ-Linolenic Acid (C20:3 n-6)", 
               "Eicosatetraenoic Acid (C20:4 n-3)", "Arachidonic Acid (C20:4 n-6)", "Eicosapentaenoic Acid (C20:5 n-3)",
               "Heneicosylic Acid (C21)", "Erucic Acid (C22:1)", "Adrenic Acid (C22:4 n-6)", "Docosapentaenoic Acid (C22:5 n-3)",
               "Osbond Acid (C22:5 n-6)", "Docosahexaenoic Acid (C22:6)", "Nervonic Acid (C24:1)", "C4", "C6", "Fat Mass",
               "Waist Circumference", "Hip Circumference", "Cholesterol", "Triglycerides", "HDL", "LDL", "VLDL", "Systolic Blood
               Pressure (SBP)", "Diastolic Blood Pressure (DBP)", "Alanine Transaminase (ALT)", "Aspartate Aminotransferase (AST)",
               "Fasting Insulin", "Fasting Glucose", "2 Hr Glucose Test", "Matsuda Index", "HOMA-IR", "Dio", "Adiponectin",
               "Height", "Weight")
snpNames = names(d)[c(49, 60, 59, 57, 47, 55, 54, 53, 52, 51, 50, 65, 48, 66, 64, 61, 63, 67, 68, 69)]

for (j in 1:length(phenoList)){
  df.temp <- FP %>% filter(Phenotype == phenoList[j])
  cc1 <- ifelse(df.temp$SNP %in% c("genors2397142", "genors9357760", "genors7744440"),
                "#5DA899", ifelse(df.temp$SNP %in% c("genors953413", "genors3798719", "genors3798713", "genors3734398",
                                                     "genors1570069", "genors2281591"), "#94CBEC", "black"))
  if (phenoList[j] == phenoEstList[j]) {
    Beta <- "%"
  } else if (phenoList[j] == phenoEstList2){
    Beta <- "%"
  } else if (phenoEstList[j] == "NA") {
    Beta <- phenoUnitList[j]
  }
  FPL <- ggplot(df.temp, aes(x=SNP, y=as.numeric(Est))) +
    geom_hline(yintercept = 0, color = "#0072b2", linewidth = 2) +
    geom_errorbar(aes(ymin=as.numeric(Conf.low), ymax=as.numeric(Conf.high)), 
                  width = 0.5, size  = 1,
                  position = "dodge", color=cc1) + coord_flip() +
    scale_x_discrete(limits=c("genors174455", "genors174602", "genors174594", "genors174576",
                              "genors1535", "genors174554", "genors174548", "genors174547", "genors174546",
                              "genors174538", "genors174537", "genors3798719", "genors1570069", "genors953413", "genors3798713",
                              "genors2281591", "genors3734398", "genors2397142", "genors9357760", "genors7744440"),
                     labels=snpNames) + 
    geom_point(size=3, color=cc1) + theme(axis.text.y = element_text(color = cc1, size = 17)) +
    theme(axis.text.x = element_text(size = 16)) + theme(axis.title.y = element_text(size = 22)) +
    labs(y=paste("Linear Model Estimate (", Beta,")")) + theme(axis.title.x = element_text(size = 20))
  print(j)
  ggsave(device="png", filename= paste("FP",phenoList[j]), path = "FPs F", width=7, height=6, units = "in")
}

#ggsave above automatically stores FP images to working directory, path defines another folder if desired

##############################################################################################################################

#Forest plots for PUFA ratios
#fatty acids should ***NOT*** be log transformed for this analysis

snpList = names(d)[c(108, 110:111, 113:114, 116:130)]
phenoList = names(d)[c(45, 46, 70:77)]

TFP<-NULL; tmp.vec<-NULL
for (i in 1:length(snpList)) {
  for (j in 1:length(phenoList)) {
    tmp.formula <- paste0(phenoList[j], " ~ age + sex + bmi + (1|ego_id) + t2d +", snpList[i])
    tmp.lm <- lmekin(as.formula(tmp.formula), data=d)
    tmp.p <- as.numeric(2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm))))))[7]
    SE <- as.numeric(sqrt(diag(vcov(tmp.lm)))[7])
    conf.low <- tmp.lm$coefficients$fixed[7] - (1.96 * SE)
    conf.high <- tmp.lm$coefficients$fixed[7] + (1.96 * SE)
    tmp.vec <- c(snpList[i], phenoList[j], as.numeric(tmp.lm$coefficients$fixed[7]), "Ratio" , conf.low, conf.high, tmp.p)
    print(snpList[i])
    print(phenoList[j])
    TFP <- rbind(TFP, tmp.vec)
  }
}

TFP <- as.data.frame(TFP)
names(TFP) = c("SNP","Phenotype","Est", "Est unit", "Conf.low", "Conf.high", "p_val")
rownames(TFP)<-NULL
View(TFP)

write.csv(TFP,"~/Chilton Laboratory/Gene-Diet Cohorts/AIR/byT2D/FP ratios.csv", row.names = TRUE)

#########################################################################################################
#loop for PUFA ratio FP

phenoList = names(d)[c(45, 46, 70:77)]
phenoNames = c("ARA/EPA Ratio", "ARA/DGLA Ratio", "ARA/DHA Ratio", "ARA/DPA n-3 Ratio", "ARA/DPA n-6 Ratio",
               "EPA/DPA n-3 Ratio", "EPA/DPA n-6 Ratio", "EPA/DHA Ratio", "DHA_DPA n-3 Ratio",
               "DHA/EPA Ratio")
snpNames = names(d)[c(49, 60, 59, 57, 47, 55, 54, 53, 52, 51, 50, 65, 48, 66, 64, 61, 63, 67, 68, 69)]

for (j in 1:length(phenoList)){
  df.temp <- TFP %>% filter(Phenotype == phenoList[j])
  cc1 <- ifelse(df.temp$SNP %in% c("genors2397142", "genors9357760", "genors7744440"),
                "#5DA899", ifelse(df.temp$SNP %in% c("genors953413", "genors3798719", "genors3798713", "genors3734398",
                                                     "genors1570069", "genors2281591"), "#94CBEC", "black"))
  FPL <- ggplot(df.temp, aes(x=SNP, y=as.numeric(Est))) +
    geom_hline(yintercept = 0, color = "#0072b2", size = 2) +
    geom_errorbar(aes(ymin=as.numeric(Conf.low), ymax=as.numeric(Conf.high)), 
                  width = 0.5,size  = 1,
                  position = "dodge", color=cc1) + coord_flip() +
    scale_x_discrete(limits=c("genors174455", "genors174602", "genors174594", "genors174576",
                              "genors1535", "genors174554", "genors174548", "genors174547", "genors174546",
                              "genors174538", "genors174537", "genors3798719", "genors1570069", "genors953413", "genors3798713",
                              "genors2281591", "genors3734398", "genors2397142", "genors9357760", "genors7744440"),
                     labels=snpNames) + 
    geom_point(size=3, color=cc1) + theme(axis.text.y = element_text(color = cc1, size = 17)) +
    theme(axis.text.x = element_text(size = 16)) + 
    theme(axis.title.y = element_text(size = 22)) +
    labs(y=paste("Linear Model Estimate (", "Ratio",")")) + 
    theme(axis.title.x = element_text(size = 20))
  print(j)
  ggsave(device="png", filename= paste("FP",phenoList[j]), path = "FPs conc", width=7, height=6, units = "in")
}


#####################################################################################################
##Spreadsheet with FA ratios FACTOR model (het/homo alleles) for final manuscript 8/24/2023

snpList = names(d)[c(47:55, 57, 59:61, 63:69)]
phenoList = names(d)[c(45, 46, 70:77)]

T2FP<-NULL; tmp.vec<-NULL
for (i in 1:length(snpList)) {
  for (j in 1:length(phenoList)) {
    tmp.formula <- paste0(phenoList[j], " ~ age + sex + bmi + (1|ego_id) + t2d + ", snpList[i])
    tmp.lm <- lmekin(as.formula(tmp.formula), data=d)
    estHET <- ((exp(tmp.lm$coefficients$fixed[7]))-1)*100
    estHOMO <- ((exp(tmp.lm$coefficients$fixed[8])-1))*100
    tmp.p <- 2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))))
    tmp.vec <- c(snpList[i], phenoList[j], as.numeric(tmp.lm$coefficients$fixed[7]), as.numeric(tmp.lm$coefficients$fixed[8]), as.numeric(tmp.p[7]), as.numeric(tmp.p[8]))
    print(snpList[i])
    print(phenoList[j])
    T2FP <- rbind(T2FP, tmp.vec)
  }
}

T2FP <- as.data.frame(T2FP)
names(T2FP) = c("SNP","Phenotype","Het Est", "RefHomo Est", "Het p val", "RefHomo p val")
rownames(T2FP)<-NULL
View(T2FP)
