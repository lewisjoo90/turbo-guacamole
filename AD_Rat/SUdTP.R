data <- read.csv("~/Desktop/Manuscripts in Preparation/Adrienne_SU_data/mural-ET1-TT-analysis-April15.csv")

# PaCO2 — Partial pressure of carbon dioxide at sea level (765 mmHg) in arterial blood is between 35 mmHg and 45 mmHg.[9]
# PvCO2 — Partial pressure of carbon dioxide at sea level in venous blood is between 40 mmHg and 50 mmHg.[9]

# explain SU effects in WT by the (low level) pericyte activation due to the cranial window implanatation; reference Gan's work shoing the inflammatory response (eg GFAP) due to craniotomy

# explain hypercapnia induced TT elongation as indicative of individual vessel's frop in flow with rising CO2; claim that this happens under physiological conds, too (and is thus in WTs as well), 
# but that the net vascular response is that of flow increase (and correpsondingly, TT shortening) with CO2 rise (even though in a few subjects we see only positive CO2 changes though)

# R HINT:
#x[x$a>thresh]<-NA
#mean(x,na.rm=TRUE)

library(nlme)
library(ggplot2)
library(outliers)
library(lme4)

#consider these categorical variables
data$subject_id <- factor(data$subject_ID)
data$genotype <- factor(data$genotype,levels=c("wt","tg"))
data$stroke <- factor(data$stroke,levels=c("PBS","ET1"))
data$tx <- factor(data$tx,levels=c("DMSO","SU6668"))
data$ventilated <- factor(data$ventilated,levels=c("N","Y"))
data$vessel_number <- factor(data$vessel_number)
data$bolus_number <- factor(data$bolus_number)
# B = branch; M = maybe; TBD,ExitsXYZ,NotinXYZ = unknown
data$vessel_type <-factor(data$vessel_type,levels=c("A","AB","Cap","CapM","V","VB","AM","VM","ABM","VBM"))

# consider ONLY ventilated animals (the TT data from non-ventilated ones appears to be pure noise, tcm4 tCO2 measurements notwithstanding!)
nrow(data) #1851
data <-subset(data,data$ventilated=="Y")
nrow(data) #1236

# looked at nonvent_data only for the purpose of doing penetrating vessel segmentation in some of them as replacements for missing XYZs of ventilated animals
# nonvent_data <-subset(data,data$ventilated=="N")
# length(unique(nonvent_data$subject_id))
# #18
# nonvent_data <-subset(nonvent_data,nonvent_data$quality!=3)
# length(unique(nonvent_data$subject_id))
# #15
# unique(nonvent_data$subject_id)
# for (i in (unique(nonvent_data$subject_id))) {
#   print(i)
#   print(unique(nonvent_data[nonvent_data$subject_id==i,]$genotype))
#   print(unique(nonvent_data[nonvent_data$subject_id==i,]$tx))
# }
# 
# nonvent_data <-subset(nonvent_data,nonvent_data$quality==3)
# length(unique(nonvent_data$subject_id))
# #3
# unique(nonvent_data$subject_id)
# for (i in (unique(nonvent_data$subject_id))) {
#   print(i)
#   print(unique(nonvent_data[nonvent_data$subject_id==i,]$genotype))
#   print(unique(nonvent_data[nonvent_data$subject_id==i,]$tx))
# }



length(unique(data$subject_id))
#26 (from 44 total; ie 18 non-ventilated data sets were acquired!!!)

#badbol_data <-subset(data,data$quality==3)
#print(unique(badbol_data$subject_id))

# consider only decent quality data
data <-subset(data,data$quality!=3)
nrow(data) #1091 ie lost 145 vessels due to poor image quality

length(unique(data$subject_id))
#22
unique(data$subject_id)

# regroup vessel_type levels
print(data$vessel_type)
sum(is.na(data$vessel_type))
#86 (out of total 1091)

levels(data$vessel_type)<-list("Aall"=c("A","AB","AM","ABM"),"Vall"=c("V","VB","VM","VBM"),"Call"=c("Cap","CapM"))
levels(data$vessel_type)
print(data$vessel_type)


# look for OUTLIERS and remove them before centering the data!!
# outliers need to be identified and then removed within coherent cohorts, ie within: wt_dmso, wt_su6668, tg_dmso,tg_su6668
# similarly, mean needs to be removed within each one of the 4 groups separately

data$TP_rest_mad_score <- NA
data$TP_CO2_mad_score <- NA

data[data$genotype=="wt" & data$tx=="DMSO",]$TP_rest_mad_score = scores(subset(data,genotype=="wt"&tx=="DMSO")$TP_rest,type="mad")
data[data$genotype=="wt" & data$tx=="SU6668",]$TP_rest_mad_score = scores(subset(data,genotype=="wt"&tx=="SU6668")$TP_rest,type="mad")
data[data$genotype=="tg" & data$tx=="DMSO",]$TP_rest_mad_score = scores(subset(data,genotype=="tg"&tx=="DMSO")$TP_rest,type="mad")
data[data$genotype=="tg" & data$tx=="SU6668",]$TP_rest_mad_score = scores(subset(data,genotype=="tg"&tx=="SU6668")$TP_rest,type="mad")

data[data$genotype=="wt" & data$tx=="DMSO",]$TP_CO2_mad_score = scores(subset(data,genotype=="wt"&tx=="DMSO")$TP_CO2,type="mad")
data[data$genotype=="wt" & data$tx=="SU6668",]$TP_CO2_mad_score = scores(subset(data,genotype=="wt"&tx=="SU6668")$TP_CO2,type="mad")
data[data$genotype=="tg" & data$tx=="DMSO",]$TP_CO2_mad_score = scores(subset(data,genotype=="tg"&tx=="DMSO")$TP_CO2,type="mad")
data[data$genotype=="tg" & data$tx=="SU6668",]$TP_CO2_mad_score = scores(subset(data,genotype=="tg"&tx=="SU6668")$TP_CO2,type="mad")

nrow(data) # 1091
# this may be too conservative: Q: should we allow unpaired measurements so that we have more of them in the end
data = subset(data, TP_rest_mad_score<3.5 & TP_CO2_mad_score<3.5)
nrow(data) # 1008 ie 83 vessels are lost due to the above outlier detection and removal scheme

sum(!is.na(data$vessel_type)) # 925 vessels

# look at TCM4 CO2 TENSION measurements
range(subset(data,genotype=="tg")$tCO2_rest)
# 29 63
range(subset(data,genotype=="wt")$tCO2_rest)
# 26 57

# TAKEAWAY: no significant effect of either genotype or treatment on resting CO2 tension
summary(lme(tCO2_rest ~ tx*genotype, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
#txSU6668             1.04799  5.851645  18  0.179094  0.8599
#genotypetg           4.19805  6.140713  18  0.683642  0.5029

# TAKEAWAY: no significant effect of either genotype or treatment on maximum CO2 tension reached during the hypercapnic challenge
summary(lme(tCO2max_CO2 ~ tx*genotype, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
#txSU6668            -10.84293  6.325920  18 -1.714047  0.1037
#genotypetg            8.17978  6.645903  18  1.230800  0.2342

# CENTER DATA, removing the per-subject-per-bolus mean TP in each of the 4 groups for each of TP_rest, TP_CO2
data$TP_rest_centered <- NA
data$TP_CO2_centered <- NA

 for (i in levels(data$subject_id)) {
   for (k in unique(data[data$subject_id==i,]$bolus_number)) {
     mrest = mean(data[data$subject_id==i & data$bolus_number==k,]$TP_rest,na.rm=TRUE)
     mCO2 = mean(data[data$subject_id==i & data$bolus_number==k,]$TP_CO2,na.rm=TRUE)
     data[data$subject_id==i & data$bolus_number==k,]$TP_rest_centered = data[data$subject_id==i & data$bolus_number==k,]$TP_rest - mrest
     data[data$subject_id==i & data$bolus_number==k,]$TP_CO2_centered = data[data$subject_id==i & data$bolus_number==k,]$TP_CO2 - mCO2
     #print(diff(range(data[data$subject_id==i & data$bolus_number==k,]$TP_rest)))
     #print(diff(range(data[data$subject_id==i & data$bolus_number==k,]$TP_rest_centered)))      
   }
 } 
 
summary(data$TP_rest_centered)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-17.6000  -1.3160   0.1047   0.0000   1.1670  20.9100 
summary(data$TP_CO2_centered)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-7.8670 -1.5020 -0.1264  0.0000  1.3460 12.4800 

# VARIANCE TESTS, COEFFICIENT OF VARIATION CONSIDERATIONS - confusing so abolished (possibly due to variability in resting levels of CO2 tension - so confounded; safer to look at HC effects)
# look at the effect of genotype, treatment on dispersion in TP_rest_centered

var.test(data[data$genotype=="wt" & data$tx=="SU6668",]$TP_rest_centered,data[data$genotype=="wt" & data$tx=="DMSO",]$TP_rest_centered)
# F test to compare two variances
# 
# data:  data[data$genotype == "wt" & data$tx == "SU6668", ]$TP_rest_centered and data[data$genotype == "wt" & data$tx == "DMSO", ]$TP_rest_centered 
# F = 0.0913, num df = 287, denom df = 257, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1 
# 95 percent confidence interval:
#   0.07186689 0.11580478 
# sample estimates:
#   ratio of variances 
# 0.09130044 smaller dispersion in wt su than wt dmso - confusing since former more reactive
x11(); plot(data[data$genotype=="wt",]$tx,data[data$genotype=="wt",]$TP_rest_centered,xlab='Treatment',ylab='Transit Time')

var.test(data[data$genotype=="tg" & data$tx=="DMSO",]$TP_rest_centered,data[data$genotype=="wt" & data$tx=="DMSO",]$TP_rest_centered)
# 
# F test to compare two variances
# 
# data:  data[data$genotype == "tg" & data$tx == "DMSO", ]$TP_rest_centered and data[data$genotype == "wt" & data$tx == "DMSO", ]$TP_rest_centered 
# F = 0.0534, num df = 242, denom df = 257, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1 
# 95 percent confidence interval:
#   0.04166004 0.06852766 
# sample estimates:
#   ratio of variances 
# 0.05340575 smaller dispersion of TP_rest_centered in tg dmso vs wt dmso - confusing since former less reactive
x11(); plot(data[data$tx=="DMSO",]$genotype,data[data$tx=="DMSO",]$TP_rest_centered,xlab='Genotype',ylab='Transit Time')


var.test(data[data$genotype=="tg" & data$tx=="SU6668",]$TP_rest_centered,data[data$genotype=="tg" & data$tx=="DMSO",]$TP_rest_centered)
# F test to compare two variances
# 
# data:  data[data$genotype == "tg" & data$tx == "SU6668", ]$TP_rest_centered and data[data$genotype == "tg" & data$tx == "DMSO", ]$TP_rest_centered 
# F = 2.5647, num df = 218, denom df = 242, p-value = 1.707e-12
# alternative hypothesis: true ratio of variances is not equal to 1 
# 95 percent confidence interval:
#   1.980300 3.327394 
# sample estimates:
#   ratio of variances 
# 2.564669  # LARGER VARIANCE OF CENTERED TP_REST values in SU tgs than in DMSO tgs - confusing since the former show higher reactivity
x11(); plot(data[data$genotype=="tg",]$tx,data[data$genotype=="tg",]$TP_rest_centered,xlab='Treatment',ylab='Tg Transit Time During Normocapnia [s]')

# consider effects of genotype, tx on coefficient of variation in TP_rest (sigma/mu)
library(raster)
cv_wt_dmso = cv(data[data$genotype=="wt" & data$tx=="DMSO",]$TP_rest,na.rm=TRUE)
print(cv_wt_dmso) #47
cv_wt_su = cv(data[data$genotype=="wt" & data$tx=="SU6668",]$TP_rest,na.rm=TRUE)
print(cv_wt_su) #30 - lower coefficient of variation in su wt than in dmso wt - confusing

cv_tg_dmso = cv(data[data$genotype=="tg" & data$tx=="DMSO",]$TP_rest,na.rm=TRUE)
print(cv_tg_dmso) # 25
cv_tg_su = cv(data[data$genotype=="tg" & data$tx=="SU6668",]$TP_rest,na.rm=TRUE)
print(cv_tg_su) # 40 - higher coefficient of variation in su tg than in dmso tg


# ALTERNATIVE WAY TO CENTER data: residuals from per-subject*bolus_number lm fit of TP data are the per-subject-per-bolus centered (mean removed :) ) times to peak (which is perhaps more robust than TPi-min(TPi) normalization) 
#data[data$genotype=="wt" & data$tx=="DMSO",]$TP_rest_centered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="wt"&tx=="DMSO")))
#data[data$genotype=="wt" & data$tx=="SU6668",]$TP_rest_centered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="wt"&tx=="SU6668")))
#data[data$genotype=="tg" & data$tx=="DMSO",]$TP_rest_centered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO")))
#data[data$genotype=="tg" & data$tx=="SU6668",]$TP_rest_centered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668")))
#
#data[data$genotype=="wt" & data$tx=="DMSO",]$TP_CO2_centered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="wt"&tx=="DMSO")))
#data[data$genotype=="wt" & data$tx=="SU6668",]$TP_CO2_centered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="wt"&tx=="SU6668")))
#data[data$genotype=="tg" & data$tx=="DMSO",]$TP_CO2_centered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO")))
#data[data$genotype=="tg" & data$tx=="SU6668",]$TP_CO2_centered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668")))

x11(); qplot(TP_rest_centered, TP_CO2_centered, colour = tx, data=subset(data,genotype=="wt")) + xlab("nTg TTP during air") + ylab("nTg TTP during hypercapnia") + layer(geom = "point", size = 3) + layer(geom = "smooth", method = "lm", size = 2) + scale_colour_hue("Treatment")
#ggsave("WT_TP_CO2_vs_TP_air.pdf", width=4, height = 3)
dev.copy(png,'WT_TP_CO2_vs_TP_air.png') # TAKEAWAY: use this in the paper as a demonstration of the SU induced decreased in HC response in the WT vessels
dev.off()

x11(); qplot(TP_rest_centered, TP_CO2_centered, colour = tx, data=subset(data,genotype=="tg")) + xlab("Tg TTP during air") + ylab("Tg TTP during hypercapnia") + layer(geom = "point", size = 3) + layer(geom = "smooth", method = "lm", size = 2) + scale_colour_hue("Treatment")
#ggsave("TG_TP_CO2_vs_TP_air.pdf", width=4, height = 3)
dev.copy(png,'TG_TP_CO2_vs_TP_air.png')
dev.off()

# get Rsquared on these regressions (for the Figure caption)
summary(lm(TP_CO2_centered ~ TP_rest_centered, subset(data,genotype=="wt"&tx=="DMSO")))
#TP_rest_centered 4.228e-01  2.343e-02   18.05   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 2.14 on 256 degrees of freedom
#Multiple R-squared: 0.5599,  Adjusted R-squared: 0.5581 
#F-statistic: 325.6 on 1 and 256 DF,  p-value: < 2.2e-16 

summary(lm(TP_CO2_centered ~ TP_rest_centered, subset(data,genotype=="wt"&tx=="SU6668")))
#TP_rest_centered  8.379e-01  6.260e-02   13.39   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 1.826 on 286 degrees of freedom
#Multiple R-squared: 0.3852,  Adjusted R-squared: 0.383 
#F-statistic: 179.2 on 1 and 286 DF,  p-value: < 2.2e-16

summary(lm(TP_CO2_centered ~ TP_rest_centered, subset(data,genotype=="tg"&tx=="DMSO")))
#TP_rest_centered  9.374e-01  9.839e-02   9.527   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 2.015 on 241 degrees of freedom
#Multiple R-squared: 0.2736,  Adjusted R-squared: 0.2705 
#F-statistic: 90.76 on 1 and 241 DF,  p-value: < 2.2e-16

summary(lm(TP_CO2_centered ~ TP_rest_centered, subset(data,genotype=="tg"&tx=="SU6668")))
#TP_rest_centered  5.864e-01  6.956e-02   8.431 4.78e-15 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 2.165 on 217 degrees of freedom
#Multiple R-squared: 0.2467,  Adjusted R-squared: 0.2433 
#F-statistic: 71.08 on 1 and 217 DF,  p-value: 4.784e-15 

# flow change factors (and std errors on it)
# > 1/4.228e-01
# [1] 2.365184
# > 2.343e-02/(4.228e-01)^2
# > 1/8.379e-01
# [1] 1.19346
# > 6.260e-02/(8.379e-01)^2
# [1] 0.08916408
# [1] 0.1310697
# > 1/9.374e-01
# [1] 1.06678
# > 9.839e-02/(9.374e-01)^2
# [1] 0.1119698
# > 1/5.864e-01
# [1] 1.705321
# > 6.956e-02/(5.864e-01)^2
# [1] 0.2022887

# sanity check: slopes reported with lm on centered data above agree (to within the 2 sig figs) with slopes reported by lme on uncentered data below - CHECKS yippeee
summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt"&tx=="DMSO")))
#TP_rest     0.416550 0.0233125 248 17.868101       0

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt"&tx=="SU6668")))
#TP_rest     0.820000 0.0613698 277 13.361623  0.0000

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO")))
#TP_rest     0.948742 0.0980163 234 9.679432   0e+00

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&tx=="SU6668")))
# TP_rest     0.577181 0.0689776 210 8.367654  0.0000


## LME PLOTTING fun
##mod=lme(TP_CO2 ~ TP_rest, random = ~1 | bolus_number, na.action = na.omit, subset(data,genotype=="wt"&tx=="DMSO"))
##for (i in 1:nrow(coef(mod))) {
##print(i)
##x11();plot.lme(mod,TP_CO2~TP_rest|bolus_number[i],abline=c(coef(mod)[i,1],coef(mod)[i,2])) # wrinky dink as won't allow variable in bolus_number[i]
##x11();plot.lme(mod,TP_CO2~TP_rest|bolus_number[9],abline=c(coef(mod)[i,1],coef(mod)[i,2]))
##}

#slice = subset(data,genotype=="wt"&tx=="DMSO")
#mmod=lme(TP_CO2 ~ TP_rest, random = ~1 | bolus_number, na.action = na.omit, slice)
#for (i in 1:nrow(coef(mmod))) {
#  x11(); plot(slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_rest,slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_CO2)
#  abline(a=coef(mmod)[i,1],b=coef(mmod)[i,2],col="blue")
#}

# WT

# based on the LME call below: in WT, no slope differences among different vessel_type within tx
summary(lme(TP_CO2 ~ tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
#vessel_typeVall:TP_rest           0.021487 0.0810680 497  0.265050  0.7911
#vessel_typeCall:TP_rest           0.006078 0.0782219 497  0.077698  0.9381

summary(lme(TP_CO2 ~ tx*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
#TP_rest           0.418899 0.0217483 525 19.261190  0.0000
#txSU6668:TP_rest  0.392116 0.0689323 525  5.688425  0.0000 # SU6668 had a significant effect on slope in nTg

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
#txDMSO:TP_rest   0.418899 0.0217483 526 19.261190  0.0000
#txSU6668:TP_rest 0.811015 0.0654115 526 12.398648  0.0000 - SU6668 made the slope higher (ie flow factor increase lower) in nTg

summary(lme(TP_CO2 ~ vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt"&tx=="DMSO")))
#vessel_typeVall:TP_rest 0.022824 0.0841899 225 0.271107  0.7866
#vessel_typeCall:TP_rest 0.007232 0.0812537 225 0.089010  0.9292 - no effect of vessel_type in nTg DMSO

summary(lme(TP_CO2 ~ vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt"&tx=="SU6668")))
#vessel_typeVall:TP_rest -0.058092 0.1041044 272 -0.558021  0.5773
#vessel_typeCall:TP_rest -0.079990 0.0968479 272 -0.825938  0.4096 - no effect of vessel_type in nTg SU6668

# WT vs TG
summary(lme(TP_CO2 ~ genotype*tx*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
#genotypetg:TP_rest           0.533218 0.1006020 969  5.300271  0.0000 - tg have different slope than wt
#txSU6668:TP_rest             0.395265 0.0711629 969  5.554378  0.0000

summary(lme(TP_CO2 ~ genotype*tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
#genotypewt:txDMSO:TP_rest    0.419472 0.0223355 970 18.780499  0.0000 - tg dmso have higher slope (smaller HC-induced flow factor) than nTg
#genotypetg:txDMSO:TP_rest    0.952690 0.0980913 970  9.712279  0.0000
#genotypewt:txSU6668:TP_rest  0.814737 0.0675668 970 12.058237  0.0000
#genotypetg:txSU6668:TP_rest  0.573327 0.0635572 970  9.020643  0.0000

# TG

# based on the LME calls below (but not LM but LM is sensitivity starved so whatever): in Tg, should relevel A vs VC since strong trend toward vessel_type effects
summary(lme(TP_CO2 ~ tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
#vessel_typeVall:TP_rest          -0.503940 0.3301918 373 -1.526204  0.1278
#vessel_typeCall:TP_rest          -0.492655 0.2823000 373 -1.745146  0.0818
summary(lme(TP_CO2 ~ tx*vessel_type/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
#txDMSO:vessel_typeAall:TP_rest    1.251058 0.2652717 374  4.716139  0.0000
#txSU6668:vessel_typeAall:TP_rest  0.581533 0.1504679 374  3.864834  0.0001 - su increases flow response in Tg arterioles 
#txDMSO:vessel_typeVall:TP_rest    0.747118 0.2376679 374  3.143538  0.0018
#txSU6668:vessel_typeVall:TP_rest  0.376840 0.0932969 374  4.039149  0.0001 - su increases flow response in Tg venules
#txDMSO:vessel_typeCall:TP_rest    0.758403 0.1616559 374  4.691467  0.0000
#txSU6668:vessel_typeCall:TP_rest  0.443298 0.0711182 374  6.233258  0.0000 - su increases flow response in Tg capillaries

#> 1/1.251058
#[1] 0.7993235
#> 0.2652717/(1.251058)^2
#[1] 0.1694869

summary(lme(TP_CO2 ~ tx*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
#txSU6668:TP_rest -0.375804 0.1206677 444 -3.114370  0.0020 - SU affected Tgs' vessels' slope

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
#txDMSO:TP_rest   0.950917 0.1012004 445 9.396383  0.0000
#txSU6668:TP_rest 0.575114 0.0657205 445 8.750906  0.0000 - SU decreased Tgs' vessels' slope (ie increased their flow increase to HC)

# NO AGREEMENT
# based on the calls below: in WT, no slope differences among different vessel_type within tx
# summary(lm(TP_CO2_centered ~ tx*vessel_type*TP_rest_centered, na.action = na.omit, subset(data,genotype=="wt")))
#vessel_typeVall:TP_rest_centered           0.05191    0.07961   0.652 0.514646    
#vessel_typeCall:TP_rest_centered           0.11143    0.07472   1.491 0.136469 

#summary(lm(TP_CO2_centered ~ tx*vessel_type*TP_rest_centered, na.action = na.omit, subset(data,genotype=="tg")))
#vessel_typeVall:TP_rest_centered           0.21659    0.46153   0.469   0.6391  
#vessel_typeCall:TP_rest_centered          -0.16351    0.39200  -0.417   0.6768  
#summary(lm(TP_CO2_centered ~ tx*vessel_type/(1+TP_rest_centered), na.action = na.omit, subset(data,genotype=="tg")))
#txDMSO:vessel_typeAall:TP_rest_centered    0.85411    0.35014   2.439 0.015162 *  
#txSU6668:vessel_typeAall:TP_rest_centered  0.84123    0.30357   2.771 0.005856 ** 
#txDMSO:vessel_typeVall:TP_rest_centered    1.07070    0.30068   3.561 0.000416 ***
#txSU6668:vessel_typeVall:TP_rest_centered  0.42063    0.17002   2.474 0.013790 *  
#txDMSO:vessel_typeCall:TP_rest_centered    0.69060    0.17626   3.918 0.000106 ***
#txSU6668:vessel_typeCall:TP_rest_centered  0.40501    0.07924   5.111 5.04e-07 ***

# BASED on LME results, in light of the similarity of venular and capillary slopes, relevel them as one category
levels(data$vessel_type)
levels(data$vessel_type)<-list("A"=c("Aall"),"VC"=c("Vall","Call"))
levels(data$vessel_type)
print(data$vessel_type)

# RECENTER within arteriolar, venules/caps in Tgs DMSO, SU6668 during air, hypercapnia
data$TP_rest_recentered <- NA
data$TP_CO2_recentered <- NA

# must cut out rows with <NA> as vessel_type (which used to be TBD) since the RH side of the equation matching below busted otherwise
data=subset(data,!is.na(data$vessel_type))

for (i in levels(data$subject_id)) {
  for (k in unique(data[data$subject_id==i,]$bolus_number)) {
    for (j in unique(data[data$subject_id==i & data$bolus_number==k,]$vessel_type)) {
      mrest = mean(data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_rest,na.rm=TRUE)
      mCO2 = mean(data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_CO2,na.rm=TRUE)
      data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_rest_recentered = data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_rest - mrest
      data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_CO2_recentered = data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_CO2 - mCO2
      #print(diff(range(data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_rest)))
      #print(diff(range(data[data$subject_id==i & data$bolus_number==k & data$vessel_type==j,]$TP_rest_recentered)))      
    }
  } 
}

#data[data$genotype=="tg" & data$tx=="DMSO" & data$vessel_type=="A",]$TP_rest_recentered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
#data[data$genotype=="tg" & data$tx=="SU6668" & data$vessel_type=="A",]$TP_rest_recentered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="A")))
#data[data$genotype=="tg" & data$tx=="DMSO" & data$vessel_type=="VC",]$TP_rest_recentered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="VC")))
#data[data$genotype=="tg" & data$tx=="SU6668" & data$vessel_type=="VC",]$TP_rest_recentered = residuals(lm(TP_rest~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="VC")))

#data[data$genotype=="tg" & data$tx=="DMSO" & data$vessel_type=="A",]$TP_CO2_recentered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
#data[data$genotype=="tg" & data$tx=="SU6668" & data$vessel_type=="A",]$TP_CO2_recentered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="A")))
#data[data$genotype=="tg" & data$tx=="DMSO" & data$vessel_type=="VC",]$TP_CO2_recentered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="VC")))
#data[data$genotype=="tg" & data$tx=="SU6668" & data$vessel_type=="VC",]$TP_CO2_recentered = residuals(lm(TP_CO2~subject_id*bolus_number,subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="VC")))

summary(data[data$genotype=="tg",]$TP_rest_recentered)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-11.18000  -0.78580   0.07817   0.00000   0.73060   9.86300 

summary(data[data$genotype=="tg",]$TP_CO2_recentered)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-7.21100 -1.14600 -0.08877  0.00000  1.03600  8.57400 

x11(); qplot(TP_rest_recentered, TP_CO2_recentered, colour = tx, data=subset(data,genotype=="tg" & vessel_type=="A"), geom=c("point","smooth"), method="lm",xlab='Tg arteriolar TTP during air',ylab='Tg arteriolar TTP during hypercapnia')
dev.copy(png,'TG_arteriolar_TP_CO2_vs_TP_air.png')
dev.off()

x11(); qplot(TP_rest_recentered, TP_CO2_recentered, colour = tx, data=subset(data,genotype=="tg" & vessel_type=="VC"), geom=c("point","smooth"), method="lm",xlab='Tg venular/capillary TTP during air',ylab='Tg venular/capilary TTP during hypercapnia')
dev.copy(png,'TG_vencap_TP_CO2_vs_TP_air.png')
dev.off()

# A - with lm
summary(lm(TP_CO2_recentered ~ tx/(1+TP_rest_recentered)-1, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A")))
#txDMSO:TP_rest_recentered   5.859e-01  2.910e-01   2.013   0.0489 *  # in dmso TG arterioles: flow increase by a factor of 1.7
#txSU6668:TP_rest_recentered 8.229e-01  1.953e-01   4.213 9.25e-05 *** # in su6668 TG arterioles: flow increase by a factor of 1.2 - su decreases HC-induced arteriolar flow increase (consistent with more CAA on them)

# VC - with lm
summary(lm(TP_CO2_recentered ~ tx/(1+TP_rest_recentered)-1, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="VC")))
#txDMSO:TP_rest_recentered   9.232e-01  1.479e-01   6.241 1.32e-09 *** # in dmso TG venules/caps: flow increase by a factor of 1.08
#txSU6668:TP_rest_recentered 4.311e-01  7.063e-02   6.103 2.88e-09 *** # in su6668 TG venules/caps: flow increase by a factor of 2.3 - su increases HC-induced arteriolar flow increase

# A, VC - with lme
summary(lme(TP_CO2 ~ tx*vessel_type/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
#txDMSO:vessel_typeA:TP_rest     1.3056852 0.2671306 378  4.887816  0.0000 - sort of ok (since in agreement to within the size of the std error) - no change in flow
#txSU6668:vessel_typeA:TP_rest   0.6355754 0.1497646 378  4.243829  0.0000 - ok
#txDMSO:vessel_typeVC:TP_rest    0.8753178 0.1397411 378  6.263852  0.0000 - ok
#txSU6668:vessel_typeVC:TP_rest  0.4664240 0.0666205 378  7.001209  0.0000 - ok

# > 1/1.3056852
# [1] 0.7658814
# > 0.2671306/(1.3056852)^2
# [1] 0.1566919
# > 1/0.6355754
# [1] 1.573377
# > 0.1497646/(0.6355754)^2
# [1] 0.3707448
# > 1/0.8753178
# [1] 1.142442
# > 0.1397411/(0.8753178)^2
# [1] 0.1823865
# > 1/0.4664240
# [1] 2.143972
# > 0.0666205/(0.4664240)^2
# [1] 0.3062289

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A")))
#txDMSO:TP_rest   1.020210 0.2383318 43 4.280628  0.0001
#txSU6668:TP_rest 0.643276 0.1784976 43 3.603833  0.0008 - ok

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="VC")))
#txDMSO:TP_rest   0.894137 0.1496749 323 5.973858  0.0000 - ok
#txSU6668:TP_rest 0.430251 0.0695075 323 6.189992  0.0000 - ok

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A"&tx=="DMSO"))) 
#TP_rest     0.949717 0.2583949 27 3.675447  0.0010 - ?!?

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A"&tx=="SU6668"))) 
#TP_rest     0.5903601 0.1603095 15 3.682626  0.0022 - ok

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="VC"&tx=="DMSO"))) 
#TP_rest     0.893965 0.1443823 160 6.191651  0.0000 - ok

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="VC"&tx=="SU6668")))
#TP_rest     0.430665 0.0731781 162 5.885161   0e+00 - ok


### ENOUGH ###

slice = subset(data,genotype=="tg"&vessel_type=="A"&tx=="DMSO")
mod=lme(TP_CO2 ~ TP_rest, random = ~1 | bolus_number, na.action = na.omit, slice)
coef(mod)
x11();plot.lme(mmod,TP_CO2~TP_rest|bolus_number,abline=c(coef(mod)[1,1],coef(mod)[1,2]))


#for (i in 1:nrow(coef(mmod))) {
#  x11(); plot(slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_rest,slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_CO2)
#  abline(a=coef(mmod)[i,1],b=coef(mmod)[i,2],col="blue")
#}



# MODELS RESULTING IN 0.59 SLOPE
summary(lm(TP_CO2_recentered ~ TP_rest_recentered, na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
#TP_rest_recentered 5.859e-01  3.007e-01   1.948   0.0597 - this actually almost includes 1 (0.6+/-0.3) ie almost no flow change

summary(lm(TP_CO2_recentered ~ TP_rest_recentered, na.action = na.omit, subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="A")))
#TP_rest_recentered 8.229e-01  1.848e-01   4.453    2e-04 *** - this includes 1 (1.2+/-0.3) ie no flow change

summary(lm(TP_CO2_recentered ~ TP_rest_recentered, na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="VC")))
#TP_rest_recentered 9.084e-01  1.423e-01   6.385 1.67e-09 *** - 1.1+/-0.1 ie no flow change
  
summary(lm(TP_CO2_recentered ~ TP_rest_recentered, na.action = na.omit, subset(data,genotype=="tg"&tx=="SU6668"&vessel_type=="VC")))
#TP_rest_recentered 4.311e-01  7.291e-02   5.912 1.82e-08 *** - 2.3+/-0.4 - su partially restores HC-induced venular/capilary flow increase  - now that the relative compliance of A vs VC is partially restored, by su, get more flow through venules/caps



#  lmobv<-lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A"))
  # observed versus fitted values
#  x11();plot(lmobv, TP_CO2~fitted(.), abline=c(0,1))
  #x11();plot(lmobv, fitted(.)~TP_rest)

#  slice = subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")
#  mmod=lme(TP_CO2 ~ TP_rest, random = ~1+TP_rest | bolus_number, na.action = na.omit, slice)
#  for (i in 1:nrow(coef(mmod))) {
#    x11(); plot(slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_rest,slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_CO2)
#    abline(a=coef(mmod)[i,1],b=coef(mmod)[i,2],col="blue")
#  }

m= lmer(TP_CO2 ~ TP_rest + (TP_rest|bolus_number), na.action = na.omit, slice)

#  mmod=lmer(TP_CO2 ~ TP_rest + (1|bolus_number), na.action = na.omit, slice)
#  slinky=coef(mmod)
#  for (i in 1:nrow(slinky$bolus_number)) {
#    x11(); plot(slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_rest,slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_CO2)
#    abline(a=coef(mmod)[i,1],b=coef(mmod)[i,2],col="blue")
#  }


summary(lmer(TP_CO2 ~ TP_rest + (1|subject_id) + (1|bolus_number), na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
# TP_rest       0.9497     0.2584   3.675

summary(lmer(TP_CO2 ~ TP_rest + (1|subject_id) + (1|subject_id:bolus_number), na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
#TP_rest       0.9497     0.2584   3.675

# slope sort of close to lm call
#summary(lmer(TP_CO2 ~ TP_rest + (1+TP_rest|subject_id)+(1+TP_rest|bolus_number), na.action = na.omit, subset(data,genotype=="tg"&tx=="DMSO"&vessel_type=="A")))
# TP_rest       0.6503     0.4196   1.550

# qplot
#x11(); qplot(TP_rest, TP_CO2, colour = subject_id, data=subset(data,genotype=="tg" & vessel_type=="A" & tx=="DMSO"), geom=c("point","smooth"), method="lm",xlab='Tg arteriolar TTP during air',ylab='Tg arteriolar TTP during hypercapnia')
#x11(); qplot(TP_rest, TP_CO2, colour = bolus_number, data=subset(data,genotype=="tg" & vessel_type=="A" & tx=="DMSO"), geom=c("point","smooth"), method="lm",xlab='Tg arteriolar TTP during air',ylab='Tg arteriolar TTP during hypercapnia')





### HERE - the above lm's on recentered arteriolar vessels TPs don't to agree with lme below - how come???

x11(); qplot(TP_rest_recentered, TP_CO2_recentered, colour = tx, data=subset(data,genotype=="tg" & vessel_type=="A"), geom=c("point","smooth"), method="lm",xlab='Tg arteriolar TTP during air',ylab='Tg arteriolar TTP during hypercapnia')
dev.copy(png,'TG_arteriolar_TP_CO2_vs_TP_air.png')
dev.off()


x11(); qplot(TP_rest_recentered, TP_CO2_recentered, colour = tx, data=subset(data,genotype=="tg" & vessel_type=="VC"), geom=c("point","smooth"), method="lm",xlab='Tg venular/capillary TTP during air',ylab='Tg venular/capilary TTP during hypercapnia')
dev.copy(png,'TG_vencap_TP_CO2_vs_TP_air.png')
dev.off()

# lme on individual subjects slopes then show box and whisker plots of these slopes or box and 1/slopes
# note that the bolus_numbers are unique so implicitly separating dmso from su6668
mo=lm(TP_CO2 ~ vessel_type*bolus_number/(1+TP_rest)-1, na.action = na.omit, subset(data,genotype=="tg"))
cmo=coef(mo)
# vessel_typeA                        vessel_typeVC                       bolus_number36                       bolus_number37                       bolus_number38 
# 13.30699635                           4.17960299                          -9.32463850                         -10.49061613                         -10.00715495 
# bolus_number39                       bolus_number40                       bolus_number47                       bolus_number48                       bolus_number52 
# -43.25967409                         -17.33903967                          15.30684274                         -10.22115194                          -7.06608238 
# bolus_number53                       bolus_number56                       bolus_number57                       bolus_number60                       bolus_number61 
# -9.21481975                          -5.44077836                          13.98935316                         -12.53176141                         -11.94752010 
# bolus_number62                       bolus_number63         vessel_typeVC:bolus_number36         vessel_typeVC:bolus_number37         vessel_typeVC:bolus_number38 
# -10.90168085                          -6.51005429                          10.28729239                          13.91007223                          12.10733048 
# vessel_typeVC:bolus_number39         vessel_typeVC:bolus_number40         vessel_typeVC:bolus_number47         vessel_typeVC:bolus_number48         vessel_typeVC:bolus_number52 
# 39.85245806                          16.29349836                                   NA                           6.93998650                           7.54300994 
# vessel_typeVC:bolus_number53         vessel_typeVC:bolus_number56         vessel_typeVC:bolus_number57         vessel_typeVC:bolus_number60         vessel_typeVC:bolus_number61 
# 9.82968089                           6.30278534                          -5.78038562                          11.55099562                          12.79963067 
# vessel_typeVC:bolus_number62         vessel_typeVC:bolus_number63  vessel_typeA:bolus_number35:TP_rest vessel_typeVC:bolus_number35:TP_rest  vessel_typeA:bolus_number36:TP_rest 
# 7.69195952                           9.29824176                          -0.94950854                           1.06034086                           0.56451755 
# vessel_typeVC:bolus_number36:TP_rest  vessel_typeA:bolus_number37:TP_rest vessel_typeVC:bolus_number37:TP_rest  vessel_typeA:bolus_number38:TP_rest vessel_typeVC:bolus_number38:TP_rest 
# 0.23761774                           0.14662291                          -0.11244797                           1.26921847                           0.56966332 
# vessel_typeA:bolus_number39:TP_rest vessel_typeVC:bolus_number39:TP_rest  vessel_typeA:bolus_number40:TP_rest vessel_typeVC:bolus_number40:TP_rest  vessel_typeA:bolus_number47:TP_rest 
# 8.61703050                           1.95079117                           2.46704579                           1.47136182                          -0.50705187 
# vessel_typeVC:bolus_number47:TP_rest  vessel_typeA:bolus_number48:TP_rest vessel_typeVC:bolus_number48:TP_rest  vessel_typeA:bolus_number52:TP_rest vessel_typeVC:bolus_number52:TP_rest 
# NA                           0.77305094                           1.22958235                          -0.15154088                           0.45307417 
# vessel_typeA:bolus_number53:TP_rest vessel_typeVC:bolus_number53:TP_rest  vessel_typeA:bolus_number56:TP_rest vessel_typeVC:bolus_number56:TP_rest  vessel_typeA:bolus_number57:TP_rest 
# 1.01682217                           1.30959222                           0.06506998                           0.35710339                          -1.11850191 
# vessel_typeVC:bolus_number57:TP_rest  vessel_typeA:bolus_number60:TP_rest vessel_typeVC:bolus_number60:TP_rest  vessel_typeA:bolus_number61:TP_rest vessel_typeVC:bolus_number61:TP_rest 
# 0.15989124                           0.80195435                           0.51450465                           0.80028191                           0.43836899 
# vessel_typeA:bolus_number62:TP_rest vessel_typeVC:bolus_number62:TP_rest  vessel_typeA:bolus_number63:TP_rest vessel_typeVC:bolus_number63:TP_rest 
# 0.72634107                           0.91544646                           0.95010126                           0.91238625 

slopes <- read.csv("~/Desktop/Manuscripts in Preparation/Adrienne_SU_data/slopes.csv", quote="")

slopes$subject_id <- factor(slopes$subject_id)
slopes$tx <- factor(slopes$treatment,levels=c("DMSO","SU6668"))
slopes$bolus_number <- factor(slopes$bolus_number)
slopes$vessel_type <-factor(slopes$vessel_type,levels=c("A","VC"))

slopes$flow=1/slopes$slope

slopes=subset(slopes,!is.na(slopes$slope))
slopes$slope_mad_score = scores(slopes$slope,type="mad")
slopes = subset(slopes, slope_mad_score<3.5)

summary(lme(flow ~ tx+vessel_type, random = ~1 | subject_id, na.action = na.omit, slopes))
#Fixed effects: flow ~ tx + vessel_type 
#Value Std.Error DF    t-value p-value
#(Intercept)    0.3230267  1.280333 19  0.2522989  0.8035
#txSU6668       2.1236435  1.447631  8  1.4669780  0.1806 - trend toward higher flow in SU vs DMSO treated Tg vessels
#vessel_typeVC -0.2251042  1.444411 19 -0.1558450  0.8778 - non-significantly lower flow in VC vs A

summary(lme(flow ~ tx, random = ~1 | subject_id, na.action = na.omit, slopes))
#txSU6668    2.1236435  1.422185  8 1.4932258  0.1737

aslopes=subset(slopes,slopes$vessel_type=="A")
x11();plot(aslopes$tx,aslopes$slope)

vcslopes=subset(slopes,slopes$vessel_type=="VC")
x11();plot(vcslopes$tx,vcslopes$flow)


#  mmod=lm(TP_CO2 ~ bolus_number/(1+TP_rest) + (1|bolus_number), na.action = na.omit, slice)
#  slinky=coef(mmod)
#  for (i in 1:nrow(slinky$bolus_number)) {
#    x11(); plot(slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_rest,slice[slice$bolus_number==unique(slice$bolus_number)[i],]$TP_CO2)
#    abline(a=coef(mmod)[i,1],b=coef(mmod)[i,2],col="blue")
#  }


# SYNTAX HINT: genotype*tx/(1+TP_rest)-1 (-1 - "don't fit and overall intercept")
# NOTE: forget changes in TP amplitudes due to timing uncertainty in the experiment and stronger effects on global vs local flow
# instead rely on increase in dispersion of TP across different vessels as a reporter on flow change
# HENCE vessel_number no longer a random effect (as was the case in the dtt amplitude analysis)!!!

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "wt") 
# AIC      BIC    logLik
# 1986.573 2014.854 -986.2863
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0009461516
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.438521 2.326836
# 
# Fixed effects: TP_CO2 ~ tx/(1 + TP_rest) - 1 
# Value Std.Error  DF   t-value p-value
# txDMSO           7.509549 1.0472581   8  7.170676  0.0001
# txSU6668         2.673420 1.3987643   8  1.911273  0.0923
# txDMSO:TP_rest   0.365039 0.0232925 408 15.671946  0.0000 TAKE AWAY: WT DMSO animals:  TT_CO2 = 0.3*TT_air (flow increase by a factor of 1/0.365=2.7)
# txSU6668:TP_rest 0.829894 0.0777632 408 10.672063  0.0000 TAKE AWAY: WT SU6668 animals: TT_CO2 = 0.8*TT_air (flow increase by a factor of 1/0.829=1.2)
# Correlation: 
#   txDMSO txSU6668 tDMSO:
#   txSU6668          0.000                
# txDMSO:TP_rest   -0.443  0.000         
# txSU6668:TP_rest  0.000 -0.780    0.000
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.46903005 -0.44713790 -0.05318695  0.48167343  4.29487205 
# 
# Number of Observations: 424
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           15

summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "tg") 
# AIC      BIC    logLik
# 2080.144 2109.032 -1033.072
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.001001147
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.881235 2.120392
# 
# Fixed effects: TP_CO2 ~ tx/(1 + TP_rest) - 1 
# Value Std.Error  DF  t-value p-value
# txDMSO           5.061190 1.2959457   8 3.905403  0.0045
# txSU6668         3.927199 1.2507588   8 3.139853  0.0138
# txDMSO:TP_rest   0.950917 0.1012004 445 9.396383  0.0000 # TAKE AWAY: TG DMSO   TP_CO2 = 0.95*TP_air (flow increase by 1.05 times ie really no change)
# txSU6668:TP_rest 0.575114 0.0657205 445 8.750906  0.0000 # TAKE AWAY: TG SU6668 TP_CO2 = 0.575*TP_air (flow increase by 1.7 times)
# Correlation: 
#   txDMSO txSU6668 tDMSO:
#   txSU6668          0.000                
# txDMSO:TP_rest   -0.604  0.000         
# txSU6668:TP_rest  0.000 -0.567    0.000
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.99275238 -0.51764893 -0.07565215  0.45401314  5.85546150 
# 
# Number of Observations: 462
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
#10                           16 

# look at the interaction term p value first to establish if there is a difference in slopes; then look at the slope values and corresponding p values
summary(lme(TP_CO2 ~ TP_rest*genotype*tx, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
# Linear mixed-effects model fit by REML
# Data: data 
# AIC      BIC    logLik
# 4064.722 4117.276 -2021.361
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev:  0.00106031
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.674065 2.221572
# 
# Fixed effects: TP_CO2 ~ TP_rest * genotype * tx 
# Value Std.Error  DF   t-value p-value
# (Intercept)                  7.489272 1.1169475 851  6.705124  0.0000 # intercept
# TP_rest                      0.366180 0.0222824 851 16.433630  0.0000 # slope of wt
# genotypetg                  -2.463497 1.6840777  16 -1.462817  0.1629 # intercept difference tg vs wt
# txSU6668                    -4.867504 1.8105901  16 -2.688352  0.0162 # intercept difference su vs dmso
# TP_rest:genotypetg           0.587479 0.1077773 851  5.450863  0.0000 # difference in slopes btw wt and tg - TAKE AWAY: slopes in wt, tg animals are different
# TP_rest:txSU6668             0.467350 0.0785306 851  5.951189  0.0000 # difference in slope btw su and dmso - TAKE AWAY: slopes of su vs dmso animals are different
# genotypetg:txSU6668          3.801353 2.5157144  16  1.511043  0.1503 # intercept difference tg and su vs wt and dmso
# TP_rest:genotypetg:txSU6668 -0.848705 0.1481296 851 -5.729479  0.0000 # difference in slopes tg and su vs wt and dmso - TAKE AWAY: slope of tg and su vs wt and dmso are different
# Correlation: 
#   (Intr) TP_rst gntypt tSU666 TP_rs: TP_:SU g:SU66
# TP_rest                     -0.397                                          
# genotypetg                  -0.663  0.263                                   
# txSU6668                    -0.617  0.245  0.409                            
# TP_rest:genotypetg           0.082 -0.207 -0.527 -0.051                     
# TP_rest:txSU6668             0.113 -0.284 -0.075 -0.629  0.059              
# genotypetg:txSU6668          0.444 -0.176 -0.669 -0.720  0.353  0.453       
# TP_rest:genotypetg:txSU6668 -0.060  0.150  0.384  0.333 -0.728 -0.530 -0.605
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.64301422 -0.48115806 -0.06793117  0.46968877  5.58711926 
# 
# Number of Observations: 886
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 20                           31 

# equivalently (note that slope values below are the same as above; intercept values below are 0 by design, due to centering of the data):
# summary(lme(TP_CO2_centered ~ TP_rest_centered*genotype*tx, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
# Linear mixed-effects model fit by REML
# Data: data 
# AIC      BIC    logLik
# 3938.284 3990.838 -1958.142
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 3.723852e-05
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev: 7.285496e-05 2.187184
# 
# Fixed effects: TP_CO2_centered ~ TP_rest_centered * genotype * tx 
# Value  Std.Error  DF   t-value p-value
# (Intercept)                           0.0000000 0.16124135 851  0.000000       1
# TP_rest_centered                      0.3699376 0.02207819 851 16.755792       0
# genotypetg                            0.0000000 0.21374070  16  0.000000       1
# txSU6668                              0.0000000 0.21431557  16  0.000000       1
# TP_rest_centered:genotypetg           0.5674280 0.10904773 851  5.203483       0
# TP_rest_centered:txSU6668             0.4768228 0.08090086 851  5.893914       0
# genotypetg:txSU6668                   0.0000000 0.29573839  16  0.000000       1
# TP_rest_centered:genotypetg:txSU6668 -0.8277536 0.15127790 851 -5.471742       0
# Correlation: 
#   (Intr) TP_rs_ gntypt tSU666 TP_r_: TP__:S g:SU66
# TP_rest_centered                      0.000                                          
# genotypetg                           -0.754  0.000                                   
# txSU6668                             -0.752  0.000  0.568                            
# TP_rest_centered:genotypetg           0.000 -0.202  0.000  0.000                     
# TP_rest_centered:txSU6668             0.000 -0.273  0.000  0.000  0.055              
# genotypetg:txSU6668                   0.545  0.000 -0.723 -0.725  0.000  0.000       
# TP_rest_centered:genotypetg:txSU6668  0.000  0.146  0.000  0.000 -0.721 -0.535  0.000
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.73260067 -0.48077589 -0.05391284  0.47618220  5.68347809 
# 
# Number of Observations: 886
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 20                           31 

summary(lme(TP_CO2 ~ genotype*tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
# Linear mixed-effects model fit by REML
# Data: data 
# AIC      BIC    logLik
# 4064.722 4117.276 -2021.361
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.001046637
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.674065 2.221572
# 
# Fixed effects: TP_CO2 ~ genotype * tx/(1 + TP_rest) - 1 
# Value Std.Error  DF   t-value p-value
# genotypewt                   7.489272 1.1169475  16  6.705124  0.0000
# genotypetg                   5.025775 1.2603752  16  3.987523  0.0011
# txSU6668                    -4.867504 1.8105901  16 -2.688352  0.0162
# genotypetg:txSU6668          3.801353 2.5157144  16  1.511043  0.1503
# genotypewt:txDMSO:TP_rest    0.366180 0.0222824 852 16.433630  0.0000 - TAKE AWAY: all of the slopes are different and significant
# genotypetg:txDMSO:TP_rest    0.953659 0.1054487 852  9.043821  0.0000
# genotypewt:txSU6668:TP_rest  0.833531 0.0753030 852 11.069018  0.0000
# genotypetg:txSU6668:TP_rest  0.572304 0.0682341 852  8.387365  0.0000
# Correlation: 
#   gntypw gntypt tSU666 gn:SU6668 gntypw:DMSO:TP_ gntypt:DMSO:TP_ gntypw:SU6668:TP_
# genotypetg                   0.000                                                                          
# txSU6668                    -0.617  0.000                                                                   
# genotypetg:txSU6668          0.444 -0.501 -0.720                                                            
# genotypewt:txDMSO:TP_rest   -0.397  0.000  0.245 -0.176                                                     
# genotypetg:txDMSO:TP_rest    0.000 -0.646  0.000  0.324     0.000                                           
# genotypewt:txSU6668:TP_rest  0.000  0.000 -0.583  0.420     0.000           0.000                           
# genotypetg:txSU6668:TP_rest  0.000  0.000  0.000 -0.293     0.000           0.000           0.000           
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.64301422 -0.48115806 -0.06793117  0.46968877  5.58711926 
# 
# Number of Observations: 886
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 20                           31

# NEED TO LOOK AT Aall, Call, Vall effects and prey that they explain why SU reinstates Tg's response to hypercapnia!!!

# is venular accumulated
# is the slope the same in all 3 compartments
# Q: if venules more dispersed than arterioles before treatment and after tx, venular dispersion similar to arterioles; consistent with increased flow through venous compartment
# dispersion down in arteriolar side and venules preserved => perfusion increased on arteriolar side

summary(lme(TP_CO2 ~ genotype*tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
# Linear mixed-effects model fit by REML
# Data: data 
# AIC      BIC    logLik
# 3522.784 3647.883 -1734.392
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0008283955
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.495112  2.05713
# 
# Fixed effects: TP_CO2 ~ genotype * tx * vessel_type * TP_rest 
# Value Std.Error  DF   t-value p-value
# (Intercept)                                  6.679522  1.489058 733  4.485738  0.0000
# genotypetg                                  -4.182220  2.393269  16 -1.747493  0.0997
# txSU6668                                    -2.777940  2.386778  16 -1.163887  0.2615
# vessel_typeVall                              4.700520  1.654404 733  2.841217  0.0046
# vessel_typeCall                              1.864656  1.316246 733  1.416647  0.1570
# TP_rest                                      0.323419  0.088958 733  3.635614  0.0003
# genotypetg:txSU6668                          3.105334  3.428822  16  0.905656  0.3786
# genotypetg:vessel_typeVall                  -0.071965  2.897942 733 -0.024833  0.9802
# genotypetg:vessel_typeCall                   1.858933  2.195378 733  0.846748  0.3974
# txSU6668:vessel_typeVall                    -3.487421  2.357549 733 -1.479257  0.1395
# txSU6668:vessel_typeCall                    -0.873628  1.896275 733 -0.460707  0.6451
# genotypetg:TP_rest                           0.951672  0.277907 733  3.424430  0.0007 ###
# txSU6668:TP_rest                             0.376674  0.165473 733  2.276349  0.0231 ###
# vessel_typeVall:TP_rest                     -0.043314  0.096384 733 -0.449388  0.6533
# vessel_typeCall:TP_rest                     -0.021158  0.090882 733 -0.232810  0.8160
# genotypetg:txSU6668:vessel_typeVall          2.988467  3.658916 733  0.816763  0.4143
# genotypetg:txSU6668:vessel_typeCall         -0.318603  2.863244 733 -0.111273  0.9114
# genotypetg:txSU6668:TP_rest                 -1.074539  0.349521 733 -3.074317  0.0022 ###
# genotypetg:vessel_typeVall:TP_rest          -0.498371  0.353845 733 -1.408444  0.1594 ## pull venules and caps into one category
# genotypetg:vessel_typeCall:TP_rest          -0.522982  0.301089 733 -1.736971  0.0828 ##
# txSU6668:vessel_typeVall:TP_rest             0.055053  0.156080 733  0.352721  0.7244
# txSU6668:vessel_typeCall:TP_rest            -0.020068  0.144314 733 -0.139059  0.8894
# genotypetg:txSU6668:vessel_typeVall:TP_rest  0.285899  0.406298 733  0.703669  0.4819
# genotypetg:txSU6668:vessel_typeCall:TP_rest  0.429349  0.351395 733  1.221844  0.2222
# Correlation: 
#   (Intr) gntypt txSU6668 vssl_V vssl_C TP_rst gn:SU6668 gnt:_V gnt:_C txSU6668:_V txSU6668:_C gn:TP_ tSU6668:T v_V:TP v_C:TP gn:SU6668:_V gn:SU6668:_C g:SU6668:T g:_V:T
# genotypetg                                  -0.622                                                                                                                                                                
# txSU6668                                    -0.624  0.388                                                                                                                                                         
# vessel_typeVall                             -0.504  0.313  0.314                                                                                                                                                  
# vessel_typeCall                             -0.639  0.398  0.399    0.635                                                                                                                                         
# TP_rest                                     -0.731  0.455  0.456    0.603  0.772                                                                                                                                  
# genotypetg:txSU6668                          0.434 -0.698 -0.696   -0.219 -0.278 -0.317                                                                                                                           
# genotypetg:vessel_typeVall                   0.287 -0.460 -0.179   -0.571 -0.363 -0.344  0.321                                                                                                                    
# genotypetg:vessel_typeCall                   0.383 -0.624 -0.239   -0.381 -0.600 -0.463  0.436     0.536                                                                                                          
# txSU6668:vessel_typeVall                     0.353 -0.220 -0.425   -0.702 -0.446 -0.423  0.296     0.401  0.267                                                                                                   
# txSU6668:vessel_typeCall                     0.444 -0.276 -0.596   -0.441 -0.694 -0.536  0.415     0.252  0.416  0.639                                                                                            
# genotypetg:TP_rest                           0.234 -0.785 -0.146   -0.193 -0.247 -0.320  0.548     0.477  0.652  0.135       0.172                                                                                
# txSU6668:TP_rest                             0.393 -0.244 -0.814   -0.324 -0.415 -0.538  0.567     0.185  0.249  0.467       0.665       0.172                                                                    
# vessel_typeVall:TP_rest                      0.650 -0.404 -0.405   -0.874 -0.749 -0.885  0.282     0.499  0.449  0.613       0.520       0.283  0.476                                                             
# vessel_typeCall:TP_rest                      0.678 -0.422 -0.423   -0.629 -0.919 -0.925  0.295     0.359  0.551  0.441       0.638       0.296  0.498     0.864                                                   
# genotypetg:txSU6668:vessel_typeVall         -0.228  0.364  0.274    0.452  0.287  0.273 -0.444    -0.792 -0.424 -0.644      -0.412      -0.378 -0.301    -0.395 -0.284                                            
# genotypetg:txSU6668:vessel_typeCall         -0.294  0.479  0.395    0.292  0.460  0.355 -0.613    -0.411 -0.767 -0.423      -0.662      -0.500 -0.440    -0.344 -0.423  0.584                                     
# genotypetg:txSU6668:TP_rest                 -0.186  0.624  0.386    0.153  0.197  0.255 -0.791    -0.379 -0.518 -0.221      -0.315      -0.795 -0.473    -0.225 -0.236  0.480        0.660                        
# genotypetg:vessel_typeVall:TP_rest          -0.177  0.552  0.110    0.238  0.204  0.241 -0.385    -0.888 -0.557 -0.167      -0.142      -0.693 -0.130    -0.272 -0.235  0.703        0.427        0.551           
# genotypetg:vessel_typeCall:TP_rest          -0.205  0.651  0.128    0.190  0.277  0.279 -0.454    -0.494 -0.895 -0.133      -0.193      -0.819 -0.150    -0.261 -0.302  0.391        0.686        0.651      0.658
# txSU6668:vessel_typeVall:TP_rest            -0.401  0.250  0.585    0.540  0.463  0.546 -0.407    -0.308 -0.277 -0.901      -0.727      -0.175 -0.719    -0.618 -0.533  0.580        0.481        0.340      0.168
# txSU6668:vessel_typeCall:TP_rest            -0.427  0.266  0.659    0.396  0.579  0.583 -0.459    -0.226 -0.347 -0.620      -0.933      -0.187 -0.810    -0.544 -0.630  0.400        0.618        0.383      0.148
# genotypetg:txSU6668:vessel_typeVall:TP_rest  0.154 -0.481 -0.225   -0.207 -0.178 -0.210  0.537     0.773  0.485  0.346       0.279       0.604  0.276     0.237  0.205 -0.882       -0.595       -0.700     -0.871
# genotypetg:txSU6668:vessel_typeCall:TP_rest  0.175 -0.558 -0.271   -0.163 -0.238 -0.239  0.630     0.423  0.766  0.255       0.383       0.701  0.333     0.223  0.259 -0.532       -0.891       -0.817     -0.564
# g:_C:T tSU6668:_V: tSU6668:_C: g:SU6668:_V:
#   genotypetg                                                                             
# txSU6668                                                                               
# vessel_typeVall                                                                        
# vessel_typeCall                                                                        
# TP_rest                                                                                
# genotypetg:txSU6668                                                                    
# genotypetg:vessel_typeVall                                                             
# genotypetg:vessel_typeCall                                                             
# txSU6668:vessel_typeVall                                                               
# txSU6668:vessel_typeCall                                                               
# genotypetg:TP_rest                                                                     
# txSU6668:TP_rest                                                                       
# vessel_typeVall:TP_rest                                                                
# vessel_typeCall:TP_rest                                                                
# genotypetg:txSU6668:vessel_typeVall                                                    
# genotypetg:txSU6668:vessel_typeCall                                                    
# genotypetg:txSU6668:TP_rest                                                            
# genotypetg:vessel_typeVall:TP_rest                                                     
# genotypetg:vessel_typeCall:TP_rest                                                     
# txSU6668:vessel_typeVall:TP_rest             0.161                                     
# txSU6668:vessel_typeCall:TP_rest             0.190  0.808                              
# genotypetg:txSU6668:vessel_typeVall:TP_rest -0.573 -0.384      -0.311                  
# genotypetg:txSU6668:vessel_typeCall:TP_rest -0.857 -0.332      -0.411       0.695      
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.40561686 -0.54275621 -0.05857954  0.49182404  4.33570495 
# 
# Number of Observations: 784
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 20                           31 


summary(lme(TP_CO2 ~ tx*vessel_type/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "tg") 
# AIC      BIC    logLik
# 1738.206 1797.543 -854.1029
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.001206994
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.739377 1.922123
# 
# Fixed effects: TP_CO2 ~ tx * vessel_type/(1 + TP_rest) - 1 
# Value Std.Error  DF   t-value p-value
# txDMSO                            2.807376 1.9200032   8  1.462173  0.1818
# txSU6668                          2.779625 1.5828152   8  1.756127  0.1171
# vessel_typeVall                   4.362884 2.2581656 373  1.932048  0.0541
# vessel_typeCall                   3.457433 1.6915018 373  2.044002  0.0417
# txSU6668:vessel_typeVall         -0.192180 2.6445319 373 -0.072671  0.9421
# txSU6668:vessel_typeCall         -0.886978 2.0457932 373 -0.433562  0.6649 # prelim TAKE AWAYS (don't quote these, first make VC a single category!):
# txDMSO:vessel_typeAall:TP_rest    1.224920 0.2650345 373  4.621736  0.0000 # TG DMSO arterioles: TP_CO2=1.2TP_air
# txSU6668:vessel_typeAall:TP_rest  0.581554 0.1500516 373  3.875698  0.0001 # TG SU arterioles: TP_CO2=0.58TP_air
# txDMSO:vessel_typeVall:TP_rest    0.735307 0.2371000 373  3.101252  0.0021 # TG DMSO venules: TP_CO2=0.73TP_air
# txSU6668:vessel_typeVall:TP_rest  0.376842 0.0930393 373  4.050350  0.0001 # TG SU venules: TP_CO2=0.38TP_air
# txDMSO:vessel_typeCall:TP_rest    0.733512 0.1618170 373  4.532976  0.0000 # TG DMSO caps: TP_CO2=0.73TP_air
# txSU6668:vessel_typeCall:TP_rest  0.443302 0.0709239 373  6.250396  0.0000 # TG SU caps: TP_CO2=0.44TP_air
# Correlation: 
#   txDMSO txSU6668 vssl_V vssl_C txSU6668:_V txSU6668:_C tDMSO:_A tSU6668:_A tDMSO:_V tSU6668:_V: tDMSO:_C
# txSU6668                          0.000                                                                                                 
# vessel_typeVall                  -0.446  0.000                                                                                          
# vessel_typeCall                  -0.619  0.000    0.504                                                                                 
# txSU6668:vessel_typeVall          0.381 -0.256   -0.854 -0.431                                                                          
# txSU6668:vessel_typeCall          0.512 -0.350   -0.417 -0.827  0.555                                                                   
# txDMSO:vessel_typeAall:TP_rest   -0.847  0.000    0.495  0.687 -0.422      -0.568                                                       
# txSU6668:vessel_typeAall:TP_rest  0.000 -0.749    0.000  0.000  0.300       0.412       0.000                                           
# txDMSO:vessel_typeVall:TP_rest   -0.117  0.000   -0.767 -0.054  0.655       0.045       0.143    0.000                                  
# txSU6668:vessel_typeVall:TP_rest  0.000 -0.273    0.000  0.000 -0.255       0.025       0.000    0.340      0.000                       
# txDMSO:vessel_typeCall:TP_rest   -0.162  0.000   -0.086 -0.538  0.073       0.445       0.200    0.000      0.239    0.000              
# txSU6668:vessel_typeCall:TP_rest  0.000 -0.372    0.000  0.000 -0.001      -0.088       0.000    0.462      0.000    0.511       0.000  
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.1267370 -0.5581911 -0.1020942  0.4894739  4.1348964 
# 
# Number of Observations: 398
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           16

# in light of the similarity of venular and capillary slopes, relevel them as one category
levels(data$vessel_type)
levels(data$vessel_type)<-list("A"=c("Aall"),"VC"=c("Vall","Call"))
levels(data$vessel_type)
print(data$vessel_type)

summary(lme(TP_CO2 ~ genotype*tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, data))
# Linear mixed-effects model fit by REML
# Data: data 
# AIC     BIC    logLik
# 3537.998 3626.23 -1749.999
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0008584166
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.591421 2.102487
# 
# Fixed effects: TP_CO2 ~ genotype * tx * vessel_type * TP_rest 
# Value Std.Error  DF   t-value p-value
# (Intercept)                                6.198117  1.528013 741  4.056326  0.0001
# genotypetg                                -3.988967  2.460115  16 -1.621456  0.1245
# txSU6668                                  -3.246591  2.419371  16 -1.341916  0.1984
# vessel_typeVC                              1.808900  1.308660 741  1.382253  0.1673
# TP_rest                                    0.366093  0.090432 741  4.048262  0.0001
# genotypetg:txSU6668                        3.409594  3.494613  16  0.975671  0.3437
# genotypetg:vessel_typeVC                   1.462364  2.152059 741  0.679519  0.4970
# txSU6668:vessel_typeVC                    -1.158222  1.884262 741 -0.614682  0.5390
# genotypetg:TP_rest                         0.957917  0.284948 741  3.361725  0.0008 #
# txSU6668:TP_rest                           0.414642  0.165815 741  2.500628  0.0126
# vessel_typeVC:TP_rest                     -0.010502  0.091696 741 -0.114526  0.9089
# genotypetg:txSU6668:vessel_typeVC          0.850739  2.825352 741  0.301109  0.7634
# genotypetg:txSU6668:TP_rest               -1.108467  0.355604 741 -3.117142  0.0019
# genotypetg:vessel_typeVC:TP_rest          -0.452938  0.297803 741 -1.520930  0.1287
# txSU6668:vessel_typeVC:TP_rest            -0.001103  0.144584 741 -0.007632  0.9939
# genotypetg:txSU6668:vessel_typeVC:TP_rest  0.299046  0.349150 741  0.856498  0.3920
# Correlation: 
#   (Intr) gntypt txSU6668 vss_VC TP_rst gn:SU6668 gn:_VC txSU6668:_VC gn:TP_ tSU6668:T v_VC:T gn:SU6668:_VC g:SU6668:T g:_VC: tSU6668:_VC:
#   genotypetg                                -0.621                                                                                                                                 
# txSU6668                                  -0.632  0.392                                                                                                                          
# vessel_typeVC                             -0.662  0.411  0.418                                                                                                                   
# TP_rest                                   -0.724  0.450  0.457    0.807                                                                                                          
# genotypetg:txSU6668                        0.437 -0.704 -0.692   -0.289 -0.317                                                                                                   
# genotypetg:vessel_typeVC                   0.402 -0.658 -0.254   -0.608 -0.491  0.463                                                                                            
# txSU6668:vessel_typeVC                     0.460 -0.285 -0.615   -0.695 -0.561  0.426     0.422                                                                                  
# genotypetg:TP_rest                         0.230 -0.782 -0.145   -0.256 -0.317  0.551     0.692  0.178                                                                           
# txSU6668:TP_rest                           0.395 -0.245 -0.805   -0.440 -0.545  0.557     0.268  0.694        0.173                                                              
# vessel_typeVC:TP_rest                      0.690 -0.429 -0.436   -0.919 -0.950  0.302     0.559  0.638        0.302  0.518                                                       
# genotypetg:txSU6668:vessel_typeVC         -0.306  0.501  0.410    0.463  0.374 -0.637    -0.762 -0.667       -0.527 -0.463    -0.426                                             
# genotypetg:txSU6668:TP_rest               -0.184  0.627  0.375    0.205  0.254 -0.785    -0.554 -0.323       -0.801 -0.466    -0.242  0.690                                      
# genotypetg:vessel_typeVC:TP_rest          -0.213  0.683  0.134    0.283  0.293 -0.481    -0.891 -0.196       -0.863 -0.160    -0.308  0.678         0.691                        
# txSU6668:vessel_typeVC:TP_rest            -0.438  0.272  0.676    0.583  0.603 -0.468    -0.354 -0.932       -0.191 -0.840    -0.634  0.622         0.392      0.195             
# genotypetg:txSU6668:vessel_typeVC:TP_rest  0.181 -0.582 -0.280   -0.241 -0.250  0.656     0.760  0.386        0.736  0.348     0.263 -0.887        -0.856     -0.853 -0.414      
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.43450980 -0.52711666 -0.05769633  0.49256270  4.17134056 
# 
# Number of Observations: 784
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 20                           31 

summary(lme(TP_CO2 ~ tx*vessel_type/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "tg") 
# AIC      BIC    logLik
# 1737.92 1781.548 -857.9602
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0009556243
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.774044 1.942051
# 
# Fixed effects: TP_CO2 ~ tx * vessel_type/(1 + TP_rest) - 1 
# Value Std.Error  DF   t-value p-value
# txDMSO                          2.4428569 1.9368860   8  1.261229  0.2428
# txSU6668                        2.3190893 1.5860329   8  1.462195  0.1818
# vessel_typeVC                   3.0635899 1.6276214 377  1.882250  0.0606
# txSU6668:vessel_typeVC         -0.0625850 1.9852789 377 -0.031525  0.9749 # TAKEAWAYs TG
# txDMSO:vessel_typeA:TP_rest     1.2862874 0.2669202 377  4.818996  0.0000 # TG DMSO arteriole: TP_CO2=1.28TP_air ie flow decrease by a factor of 0.8
# txSU6668:vessel_typeA:TP_rest   0.6355613 0.1495041 377  4.251130  0.0000 # TG SU arteriole: TP_CO2=0.63TP_air   ie flow increase by a factor of 1.6 - SU changed Tg arteriolar flow resp from decrease to increase!!
# txDMSO:vessel_typeVC:TP_rest    0.8634776 0.1397100 377  6.180501  0.0000 # TG DMSO ven/cap: TP_CO2=0.86TP_air   ie flow increase by a factor of 1.2
# txSU6668:vessel_typeVC:TP_rest  0.4664194 0.0665034 377  7.013467  0.0000 # TG SU ven/cap: TP_CO2=0.47TP_air     ie flow increase by a factor of 2.1 - SU increased Tg ven/cap flow increase
# Correlation: 
#   txDMSO txSU6668 vss_VC txSU6668:_VC tDMSO:_A tSU6668:_A tDMSO:_V
# txSU6668                        0.000                                                          
# vessel_typeVC                  -0.657  0.000                                                   
# txSU6668:vessel_typeVC          0.539 -0.360   -0.820                                          
# txDMSO:vessel_typeA:TP_rest    -0.846  0.000    0.732 -0.600                                   
# txSU6668:vessel_typeA:TP_rest   0.000 -0.743    0.000  0.427        0.000                      
# txDMSO:vessel_typeVC:TP_rest   -0.166  0.000   -0.477  0.391        0.204    0.000             
# txSU6668:vessel_typeVC:TP_rest  0.000 -0.365    0.000 -0.082        0.000    0.456      0.000  
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.23376828 -0.54555050 -0.09497923  0.51539425  4.33304732 
# 
# Number of Observations: 398
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           16 

summary(lme(TP_CO2 ~ TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A"&tx=="DMSO")))


summary(lme(TP_CO2 ~ tx/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg"&vessel_type=="A")))


summary(lme(TP_CO2 ~ tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="tg")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "tg") 
# AIC      BIC    logLik
# 1737.92 1781.548 -857.9602
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0009414242
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.774044 1.942051
# 
# Fixed effects: TP_CO2 ~ tx * vessel_type * TP_rest 
# Value Std.Error  DF   t-value p-value
# (Intercept)                     2.4428569 1.9368860 376  1.261229  0.2080
# txSU6668                       -0.1237675 2.5034033   8 -0.049440  0.9618
# vessel_typeVC                   3.0635899 1.6276214 376  1.882250  0.0606
# TP_rest                         1.2862874 0.2669202 376  4.818996  0.0000
# txSU6668:vessel_typeVC         -0.0625850 1.9852789 376 -0.031525  0.9749
# txSU6668:TP_rest               -0.6507261 0.3059377 376 -2.126989  0.0341
# vessel_typeVC:TP_rest          -0.4228098 0.2748033 376 -1.538591  0.1247 #
# txSU6668:vessel_typeVC:TP_rest  0.2536678 0.3053234 376  0.830817  0.4066
# Correlation: 
#   (Intr) txSU6668 vss_VC TP_rst txSU6668:_VC tSU6668:T v_VC:T
# txSU6668                       -0.774                                                     
# vessel_typeVC                  -0.657  0.509                                              
# TP_rest                        -0.846  0.654    0.732                                     
# txSU6668:vessel_typeVC          0.539 -0.645   -0.820 -0.600                              
# txSU6668:TP_rest                0.738 -0.801   -0.638 -0.872  0.732                       
# vessel_typeVC:TP_rest           0.737 -0.570   -0.953 -0.867  0.781        0.757          
# txSU6668:vessel_typeVC:TP_rest -0.664  0.694    0.858  0.781 -0.930       -0.872    -0.900
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.23376828 -0.54555050 -0.09497923  0.51539425  4.33304732 
# 
# Number of Observations: 398
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           16 

summary(lme(TP_CO2 ~ tx*vessel_type/(1+TP_rest)-1, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "wt") 
# AIC      BIC    logLik
# 1797.565 1840.848 -887.7823
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0006650282
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.404017 2.255679
# 
# Fixed effects: TP_CO2 ~ tx * vessel_type/(1 + TP_rest) - 1 
# Value Std.Error  DF   t-value p-value
# txDMSO                          6.222164 1.5511415   8  4.011345  0.0039
# txSU6668                        2.981616 1.9350454   8  1.540851  0.1619
# vessel_typeVC                   1.823474 1.4028049 366  1.299877  0.1945
# txSU6668:vessel_typeVC         -1.160103 2.0205693 366 -0.574147  0.5662 # TAKEAWAYS WT: there's no difference btw A vs VC in WT; with SU, both, in proportion, increase flow less upon HC (activated pericytes constrict the microvasculature)
# txDMSO:vessel_typeA:TP_rest     0.363805 0.0969210 366  3.753624  0.0002 # WT DMSO arteriole: TP_CO2=0.36TP_air - flow increase by 2.8
# txSU6668:vessel_typeA:TP_rest   0.778286 0.1476130 366  5.272478  0.0000 # WT SU arteriole: TP_CO2=0.77TP_air - flow increase by 1.3
# txDMSO:vessel_typeVC:TP_rest    0.353663 0.0306850 366 11.525596  0.0000 # WT DMSO ven/cap: TP_CO2=0.35TP_air
# txSU6668:vessel_typeVC:TP_rest  0.766293 0.0899455 366  8.519520  0.0000 # WT SU ven/cap: TP_CO2=0.76TP_air
# Correlation: 
#   txDMSO txSU6668 vss_VC txSU6668:_VC tDMSO:_A tSU6668:_A tDMSO:_V
# txSU6668                        0.000                                                          
# vessel_typeVC                  -0.699  0.000                                                   
# txSU6668:vessel_typeVC          0.485 -0.437   -0.694                                          
# txDMSO:vessel_typeA:TP_rest    -0.765  0.000    0.808 -0.561                                   
# txSU6668:vessel_typeA:TP_rest   0.000 -0.880    0.000  0.468        0.000                      
# txDMSO:vessel_typeVC:TP_rest   -0.080  0.000   -0.393  0.273        0.113    0.000             
# txSU6668:vessel_typeVC:TP_rest  0.000 -0.525    0.000 -0.141        0.000    0.584      0.000  
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.2079007 -0.5164764 -0.0046264  0.4769701  3.8968386 
# 
# Number of Observations: 386
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           15 

# NB: look at genotype*tx term (even stronger if the drug effect in opposite directions in tg vs wt)

# NOTE: argument in support of keeping WT_SU6668 group: if same effect size due to tx in both groups then not curing the disease just modifying the animal

summary(lme(TP_CO2 ~ tx*vessel_type*TP_rest, random = ~1 | subject_id/bolus_number, na.action = na.omit, subset(data,genotype=="wt")))
# Linear mixed-effects model fit by REML
# Data: subset(data, genotype == "wt") 
# AIC      BIC    logLik
# 1797.565 1840.848 -887.7823
# 
# Random effects:
#   Formula: ~1 | subject_id
# (Intercept)
# StdDev: 0.0006683731
# 
# Formula: ~1 | bolus_number %in% subject_id
# (Intercept) Residual
# StdDev:    2.404017 2.255679
# 
# Fixed effects: TP_CO2 ~ tx * vessel_type * TP_rest 
# Value Std.Error  DF   t-value p-value
# (Intercept)                     6.222164 1.5511415 365  4.011345  0.0001
# txSU6668                       -3.240548 2.4800082   8 -1.306668  0.2276
# vessel_typeVC                   1.823474 1.4028049 365  1.299877  0.1945
# TP_rest                         0.363805 0.0969210 365  3.753624  0.0002
# txSU6668:vessel_typeVC         -1.160103 2.0205693 365 -0.574147  0.5662
# txSU6668:TP_rest                0.414481 0.1765879 365  2.347166  0.0195 # significant effect of SU in WT
# vessel_typeVC:TP_rest          -0.010142 0.0983092 365 -0.103162  0.9179 # no vessel-category variation in WT
# txSU6668:vessel_typeVC:TP_rest -0.001852 0.1550327 365 -0.011943  0.9905
# Correlation: 
#   (Intr) txSU6668 vss_VC TP_rst txSU6668:_VC tSU6668:T v_VC:T
# txSU6668                       -0.625                                                     
# vessel_typeVC                  -0.699  0.437                                              
# TP_rest                        -0.765  0.478    0.808                                     
# txSU6668:vessel_typeVC          0.485 -0.644   -0.694 -0.561                              
# txSU6668:TP_rest                0.420 -0.836   -0.443 -0.549  0.699                       
# vessel_typeVC:TP_rest           0.729 -0.456   -0.919 -0.951  0.638        0.522          
# txSU6668:vessel_typeVC:TP_rest -0.462  0.705    0.583  0.603 -0.932       -0.844    -0.634
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.2079007 -0.5164764 -0.0046264  0.4769701  3.8968386 
# 
# Number of Observations: 386
# Number of Groups: 
#   subject_id bolus_number %in% subject_id 
# 10                           15 
