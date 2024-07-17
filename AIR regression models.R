
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
library(here)

setwd("~/Chilton Laboratory/Gene-Diet Cohorts/AIR/Final AIR Datasheets")
dgf<-read_excel("AIR FA geno conc final.xlsx", na = "NA")
dp <- read_excel("AIR phenotypes.xlsx", na = "NA")
dped <- read_excel ("AIR_Pedigree_forAnalysis.xlsx")
DP <- merge(dp, dped, by = "id")
d <- merge.data.frame(dgf, DP)

#set rs ids & sex as factor

#genotype as integer
d$sex <-  as.factor(d$sex)

#currently coded by ANCESTRAL allele
# ^ except rs174538 which is flipped because of 1000G haplotypes

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
d$genors1535 <-  as.integer(d$rs1535)-1
d$genors1570069 <-  as.integer(d$rs1570069)-1
d$genors174455 <-  as.integer(d$rs174455)-1
d$genors174537 <-  as.integer(d$rs174537)-1
d$genors174538 <-  as.integer(d$rs174538)-1
d$genors174546 <-  as.integer(d$rs174546)-1
d$genors174547 <-  as.integer(d$rs174547)-1
d$genors174548 <-  as.integer(d$rs174548)-1
d$genors174554 <-  as.integer(d$rs174554)-1
d$genors174570 <-  as.integer(d$rs174570)-1
d$genors174576 <-  as.integer(d$rs174576)-1
d$genors174579 <-  as.integer(d$rs174579)-1
d$genors174594 <-  as.integer(d$rs174594)-1
d$genors174602 <-  as.integer(d$rs174602)-1
d$genors2281591 <-  as.integer(d$rs2281591)-1
d$genors7115739 <-  as.integer(d$rs7115739)-1
d$genors3734398 <-  as.integer(d$rs3734398)-1
d$genors3798713 <-  as.integer(d$rs3798713)-1
d$genors3798719 <-  as.integer(d$rs3798719)-1
d$genors953413 <-  as.integer(d$rs953413)-1
d$genors2397142 <-  as.integer(d$rs2397142)-1
d$genors9357760 <-  as.integer(d$rs9357760)-1
d$genors7744440 <-  as.integer(d$rs7744440)-1


#log transform all variables for normalization, and for consistency among reported units
phenoList = names(d)[c(4, 6:33, 35:38, 41:43, 81:102)]
phenoListLOG <- phenoList[c(1:18, 19:23, 24:29, 30:36, 37:58)]
phenoLL2 <- phenoList[c(15)]
d <- d %>% mutate(across(c(all_of(phenoListLOG)), log))

#C18_2 needs an additional transformation to become normalized
d <- d %>% mutate(across(c(all_of(phenoLL2)), log))

#kinship matrix for adjusting for relatedness
ped <- pedigree(id=dped$ego_id, dadid=dped$fa_id, momid=dped$mo_id, famid=dped$fam_id, sex=dped$pedsex)
kin <- makekinship(dped$fam_id, dped$ego_id, dped$fa_id, dped$mo_id)

##################################################################################
#######~~~~~~~~~~ Data frames by T2D status ~~~~~~~~~~##############

#log transform these because they are transformed above: hba1c >= 6.5 | fast.glu >= 126 | Glu_2hr >= 200)
#hba1c >= 1.871802 | fast.glu >= 4.836282 | Glu_2hr >= 5.298317)

#add new t2d variable to original data set AFTER log transformation
d <- d %>%
  mutate(t2d = case_when(hba1c < 1.871802 & fast.glu < 4.836282 & Glu_2hr < 5.298317 ~ 'non-T2D',
                         hba1c >= 1.740466 & hba1c <= 1.856298 |  fast.glu >= 4.60517 & fast.glu <= 4.828314 |
                           Glu_2hr >= 4.941642 & Glu_2hr <= 5.293305 ~ 'pre-T2D',
                         hba1c >= 1.871802 | fast.glu >= 4.836282 | Glu_2hr >= 5.298317 ~ 'T2D',
  ))

#for not log transformed data (for use in box plots, etc)
d <- d %>%
  mutate(t2d = case_when(hba1c < 6.5 & fast.glu < 126 & Glu_2hr < 200 ~ 'non-T2D',
                         hba1c >= 5.7 & hba1c <= 6.4 |  fast.glu >= 100 & fast.glu <= 125 |
                           Glu_2hr >= 140 & Glu_2hr <= 199 ~ 'pre-T2D',
                         hba1c >= 6.5 | fast.glu >= 126 | Glu_2hr >= 200 ~ 'T2D',
  ))
###################################################################################
#linear mixed model base code
lmekin(d$C22_4n6 ~ age + sex + bmi + genors3798719 + (1|ego_id) + t2d, data=d)
LMM <- lmekin(d$homair ~ age + sex + bmi + genors7744440 + (1|ego_id) + t2d, data=d)
2*pnorm(-abs(fixef(LMM)/ sqrt(diag(vcov(LMM)))))

#facet grid base code
d %>% filter(!is.na(rs174537)) %>% 
  droplevels() %>% ggplot(aes(x=rs174537,y=(C20_4n6))) + 
  geom_boxplot() + geom_jitter(width=.2) + facet_grid(. ~rs174537) +
  theme(axis.text.x = element_blank())

#boxplot facet grid base code
ara <- d %>% filter(!is.na(rs174537)) %>%
  droplevels() %>% ggplot(aes(y=(C20_4n6), fill = rs174537)) +  geom_boxplot() +
  scale_fill_manual(values = c("#94cbec", "#5da899", "#c26a77")) +
  geom_boxplot(lwd=1) + theme(axis.text.y = element_text(size = 17),
                              axis.text.x = element_blank()) + facet_grid(. ~rs174537) +
  theme(strip.text.x = element_text(size=18)) + theme(legend.position = "none")

######################Generic FA/Phenotype Plot Code########################
#boxplot with facet grid by t2d, or sex
fa <- d %>% filter(!is.na(rs174537)) %>% filter(!is.na(t2d)) %>%
  #filter(sex=="M") %>%
  droplevels() %>% ggplot(aes(y=(EPA_DHA), fill = t2d)) +  geom_boxplot() +
  scale_fill_manual(values = c("#94cbec", "#5da899", "#c26a77")) +
  geom_boxplot(lwd=1) + theme(axis.text.y = element_text(size = 15),
                              axis.text.x = element_blank()) + facet_grid(. ~ t2d) +
  theme(strip.text.x = element_text(size=15)) + theme(legend.position = "none") +
  theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

#Clean & Labeled FA Plot
FA <- fa + labs(y = "ARA (ng/mL)") + theme(axis.title.y = element_text(size=20)) +
               theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
              theme(legend.position = "right") + theme(axis.ticks.x = element_blank())# + ylim(0,15)
      
FA    
                                                                                                                                    
##########################################################################################################

#Using lmekin (linear mixed model accounting for random effects) to adjust for AIR Cohort relatedness
#make sure to load ped & kin at beginning, for pedigree and kinship matrix data

#########################################################################################################

#MASTER TABLE heterozygous and reference homozygous alleles (factor model) - summary stats

snpList = names(d)[c(47:55, 57, 59:61, 63:69)]
phenoList = names(d)[c(4, 6:11, 13:31, 33, 35:38, 41:43, 81:98, 100:102)]
phenoEstList = c("C12", "C14", "C14_1", "C15", "C15_1", "C16", "C16_1", "C17_1", "C18", "C18_1cis", "C18_1trans",
                 "C18_1vac", "C18_2trans", "C18_2cis", "C18_3n3", "C18_3n6", "C18_4n3", "C20", "C20_1", "C20_2", "C20_3n3",
                 "C20_3n6", "C20_4n3", "C20_4n6", "C20_5n3", "C21", "C22_1", "C22_4n6", "C22_5n3", "C22_5n6", "C22_6", "C24_1", "C4", "C6",
                 "fatmass", "wc", "hip.circum", "chol", "trg", "hdl", "ldl", "vldl", "sbp", "dbp", "alt", "ast",
                 "fast.ins", "fast.glu", "Glu_2hr", "matsuda_index", "homair", "dio", "adiponectin", "height", "weight")
phenoEstList2 = ("C18_2cis")
phenoUnitList <- rep(c("ng/ml", "% (real)", "cm", "mg/dl", "mmHg", "IU/L", "uIu/ml", "mg/dl", "index units", "ug/ml", "cm", "kg"), 
                     times=c(34, 1, 2, 5, 2, 2, 1, 2, 3, 1, 1, 1))


#for calculating z value: fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))
#for calculating std error: sqrt(diag(vcov(tmp.lm)))
#for calculating p value: 2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))))

#for M/F analyses - replace data=d below with "SFP" if wanting to specify by sex.
#also need to remove sex in tmp.formula and adjust [#] references accordingly.
#SFP <- d[which(d$sex=="F"),]

results5<-NULL; tmp.vec<-NULL; tmp.results<-NULL
for (i in 1:length(snpList)) {
  for (j in 1:length(phenoList)) {
    tmp.formula <- paste0(phenoList[j], " ~ age + sex + bmi + (1|ego_id) + t2d + ", snpList[i])
    tmp.lm <- lmekin(as.formula(tmp.formula), data=d)
    estHET <- ((exp(tmp.lm$coefficients$fixed[7]))-1)*100
    estHOMO <- ((exp(tmp.lm$coefficients$fixed[8])-1))*100
    tmp.p <- 2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))))
    if (phenoList[j] == phenoEstList[j]) {
      tmp.vec <- c(snpList[i], phenoList[j], estHET, estHOMO, "%", as.numeric(tmp.p[7]), as.numeric(tmp.p[8]))
    } else if (phenoList[j] == phenoEstList2){
      estLAhet <- (exp(exp(tmp.lm$coefficients$fixed[7]))-1)*100
      estLAhomo <- (exp(exp(tmp.lm$coefficients$fixed[8]))-1)*100
      tmp.vec <- c(snpList[i], phenoList[j], estLAhet, estLAhomo, "%", as.numeric(tmp.p[7]), as.numeric(tmp.p[8])) 
    } else if (phenoEstList[j] == "NA") {
      tmp.vec <- c(snpList[i], phenoList[j], tmp.lm$coefficients$fixed[7], tmp.lm$coefficients$fixed[8], phenoUnitList[j], as.numeric(tmp.p[7]), as.numeric(tmp.p[8]))
    }
    print(snpList[i])
    print(phenoList[j])
    tmp.results <- rbind(tmp.results, tmp.vec)
    tmp.results <- as.data.frame(tmp.results)
    names(tmp.results) = c("SNP","Phenotype","Het_Est","RefHomo_Est","Est unit", "Het_pval", "RefHomo_pval")
  }
    tmp.results$p_fdr_het <- p.adjust(tmp.results$Het_pval, method="fdr")
    tmp.results$p_fdr_althomo <- p.adjust(tmp.results$AltHomo_pval, method = "fdr")
    results5 <- rbind(results5, tmp.results)
    tmp.results <- NULL
}

results5 <- as.data.frame(results5)
names(results5) = c("SNP","Phenotype","Het_Est","RefHomo_Est", "Est unit", "Het_pval", "RefHomo_pval", "Het_pval(fdr)", "RefHomo_pval(fdr)")
p.vec1 <- results5$Het_pval
p.vec2 <-results5$AltHomo_pval
global.fdr1 <- p.adjust(p.vec1, method = "fdr")
global.fdr2 <- p.adjust(p.vec2, method = "fdr")
results5$global_fdr_Het <- global.fdr1
results5$global_fdr_AltHomo <- global.fdr2
rownames(results5)<-NULL
View(results5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#MASTER TABLE (additive model) - summary stats

snpList = names(d)[c(108:116, 118, 120:122, 124:130)]
phenoList = names(d)[c(4, 6:11, 13:31, 33, 35:38, 41:43, 81:99, 100:102)]
phenoEstList = c("C12", "C14", "C14_1", "C15", "C15_1", "C16", "C16_1", "C17_1", "C18", "C18_1cis", "C18_1trans",
                 "C18_1vac", "C18_2trans", "C18_2cis", "C18_3n3", "C18_3n6", "C18_4n3", "C20", "C20_1", "C20_2", "C20_3n3",
                 "C20_3n6", "C20_4n3", "C20_4n6", "C20_5n3", "C21", "C22_1", "C22_4n6", "C22_5n3", "C22_5n6", "C22_6", "C24_1", "C4", "C6",
                 "fatmass", "wc", "hip.circum", "chol", "trg", "hdl", "ldl", "vldl", "sbp", "dbp", "alt", "ast",
                 "fast.ins", "fast.glu", "Glu_2hr", "matsuda_index", "homair", "dio", "hba1c", "adiponectin", "height", "weight")
phenoEstList2 = ("C18_2cis")
phenoUnitList <- rep(c("ng/ml", "% (real)", "cm", "mg/dl", "mmHg", "IU/L", "uIu/ml", "mg/dl", "index units", "% (real)", "ug/ml", "cm", "kg"), 
                     times=c(34, 1, 2, 5, 2, 2, 1, 2, 3, 1, 1, 1, 1))

#for M/F analyses - replace data=d below with "SFP" if wanting to specify by sex.
#also need to remove sex in tmp.formula and adjust [#] references accordingly.
#SFP <- d[which(d$sex=="M"),]

results6<-NULL; tmp.vec<-NULL; tmp.results<-NULL
for (i in 1:length(snpList)) {
  for (j in 1:length(phenoList)) {
    tmp.formula <- paste0(phenoList[j], " ~ age + sex + bmi + (1|ego_id) + t2d + ", snpList[i])
    tmp.lm <- lmekin(as.formula(tmp.formula), data=d)
    tmp.p <- 2*pnorm(-abs(fixef(tmp.lm)/ sqrt(diag(vcov(tmp.lm)))))
    est <- ((exp(tmp.lm$coefficients$fixed[7]))-1)*100
    if (phenoList[j] == phenoEstList[j]) {
      tmp.vec <- c(snpList[i], phenoList[j], est, "%", as.numeric(tmp.p[7]))
    } else if (phenoList[j] == phenoEstList2){
      estLA <- (exp(exp(tmp.lm$coefficients$fixed[7]))-1)*100
      tmp.vec <- c(snpList[i], phenoList[j], estLA, "%", as.numeric(tmp.p[7])) 
    } else if (phenoEstList[j] == "NA") {
      tmp.vec <- c(snpList[i], phenoList[j], tmp.lm$coefficients$fixed[7], phenoUnitList[j], as.numeric(tmp.p[7]))
    }
    print(snpList[i])
    print(phenoList[j])
    tmp.results <- rbind(tmp.results, tmp.vec)
    tmp.results <- as.data.frame(tmp.results)
    names(tmp.results) = c("SNP","Phenotype","Est", "Est unit", "pval")
  }
    tmp.results$p_fdr <- p.adjust(tmp.results$pval, method="fdr")
    results6 <- rbind(results6, tmp.results)
    tmp.results <- NULL
}

results6 <- as.data.frame(results6)
names(results6) = c("SNP","Phenotype","Est", "Est unit", "pval", "pval(fdr)")
rownames(results6)<-NULL
p.vec <- results6$pval
global.fdr <- p.adjust(p.vec, method = "fdr")
results6$global_fdr <- global.fdr
View(results6)

###############################################################################################################################################
write.csv(results5,"~/Chilton Laboratory/Gene-Diet Cohorts/AIR/byT2D/results5 - (factorbyT2D)FEMALE - kin - beta, p, fdr, conc.csv", row.names = TRUE)
write.csv(results6,"~/Chilton Laboratory/Gene-Diet Cohorts/AIR/byT2D/results6 - (integerbyT2D)MALE - kin - beta, p, fdr, conc.csv", row.names = TRUE)
##############################################################################################################
#conversion summaries - values should NOT be log transformed here

ggset <- d[which(d$rs174537=="G"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$wc)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])
as.numeric(summary(ggset$chol)[4])
as.numeric(summary(ggset$trg)[4])
as.numeric(summary(ggset$hdl)[4])
as.numeric(summary(ggset$ldl)[4])
as.numeric(summary(ggset$vldl)[4])
as.numeric(summary(ggset$sbp)[4])
as.numeric(summary(ggset$dbp)[4])
as.numeric(summary(ggset$fast.ins)[4])
as.numeric(summary(ggset$fast.glu)[4])
as.numeric(summary(ggset$Glu_2hr)[4])
as.numeric(summary(ggset$matsuda_index)[4])
as.numeric(summary(ggset$homair)[4])
as.numeric(summary(ggset$dio)[4])
as.numeric(summary(ggset$adiponectin)[4])
as.numeric(summary(ggset$alt)[4])
as.numeric(summary(ggset$ast)[4])

#FAs
as.numeric(summary(ggset$C20_4n6)[4])
as.numeric(summary(ggset$C20_3n6)[4])
as.numeric(summary(ggset$C20_5n3)[4])
as.numeric(summary(ggset$C22_6)[4])
as.numeric(summary(ggset$ARA_DGLA)[4])
as.numeric(summary(ggset$ARA_EPA)[4])
as.numeric(summary(ggset$ARA_DHA)[4])
as.numeric(summary(ggset$DHA_EPA)[4])

ggset <- d[which(d$rs174537=="T"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$wc)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])
as.numeric(summary(ggset$chol)[4])
as.numeric(summary(ggset$trg)[4])
as.numeric(summary(ggset$hdl)[4])
as.numeric(summary(ggset$ldl)[4])
as.numeric(summary(ggset$vldl)[4])
as.numeric(summary(ggset$sbp)[4])
as.numeric(summary(ggset$dbp)[4])
as.numeric(summary(ggset$fast.ins)[4])
as.numeric(summary(ggset$fast.glu)[4])
as.numeric(summary(ggset$Glu_2hr)[4])
as.numeric(summary(ggset$matsuda_index)[4])
as.numeric(summary(ggset$homair)[4])
as.numeric(summary(ggset$dio)[4])
as.numeric(summary(ggset$adiponectin)[4])
as.numeric(summary(ggset$alt)[4])
as.numeric(summary(ggset$ast)[4])

ggset <- d[which(d$rs174455=="A"),]
as.numeric(summary(ggset$trg)[4])
as.numeric(summary(ggset$hdl)[4])
as.numeric(summary(ggset$ldl)[4])
as.numeric(summary(ggset$vldl)[4])
as.numeric(summary(ggset$fast.ins)[4])
as.numeric(summary(ggset$fast.glu)[4])
as.numeric(summary(ggset$Glu_2hr)[4])
as.numeric(summary(ggset$homair)[4])

ggset <- d[which(d$rs174455=="G"),]
as.numeric(summary(ggset$trg)[4])
as.numeric(summary(ggset$hdl)[4])
as.numeric(summary(ggset$ldl)[4])
as.numeric(summary(ggset$vldl)[4])
as.numeric(summary(ggset$fast.ins)[4])
as.numeric(summary(ggset$fast.glu)[4])
as.numeric(summary(ggset$Glu_2hr)[4])
as.numeric(summary(ggset$homair)[4])

ggset <- d[which(d$rs174594=="A"),]
as.numeric(summary(ggset$hdl)[4])

###females only

tmpd <- d[which(d$sex=="F"),]
ggset <- tmpd[which(tmpd$rs174576=="C"),]
as.numeric(summary(ggset$hip.circum)[4])

ggset <- tmpd[which(tmpd$rs174602=="T"),]
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174455=="A"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174455=="G"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174537=="G"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174537=="T"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

###males only
tmpd <- d[which(d$sex=="M"),]
ggset <- tmpd[which(tmpd$rs174537=="G"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174537=="T"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174455=="A"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])

ggset <- tmpd[which(tmpd$rs174455=="G"),]
as.numeric(summary(ggset$fatmass)[4])
as.numeric(summary(ggset$hip.circum)[4])
as.numeric(summary(ggset$height)[4])
as.numeric(summary(ggset$weight)[4])
