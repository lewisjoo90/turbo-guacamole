library(nlme)
library(ggplot2)
library(outliers)
library(Matrix)
library(Rcpp)
library(lme4)
library(MASS)
library(data.table)
library(grid)
library(gridExtra)
library(eigenprcomp)
library(FRB)
library(MethComp)
library(deming)
setwd("~/Dropbox/Stem2/analysis/AD_Rat")
data<-as.data.frame(read.csv("bolus_AD_Rat_all_OLrej_rev4_overlap.csv",header=T))
#data_exc <- data[-c(15,16,39,40,41,42,45,46,51,52,57,58), ]

#consider these categorical variables
data$id <- factor(data$id)
# data$trt <- factor(data$trt)
data$vessel_number <- factor(data$vessel_number)
data$bolus_number <- factor(data$bolus_number)
data$state <- factor(data$state,levels=c("air","co2"))
data$depth <- factor(data$depth)

## extracting rows with matching (air&co2) data 
# break down the data by vessel type
data_v <- subset(data,type=="v")
data_a <- subset(data,type=="a")
data_c <- subset(data,type=="c")
# capillaries have row entries that are not coupled as air&co2
data_c_air <- subset(data_c, state=="air")
data_c_co2 <- subset(data_c, state=="co2")
keys <- c("vessel_number", "depth")
data_c_air_match <- merge(data_c_air, data_c_co2,by=keys)
ves <- data_c_air_match$vessel_number
data_c_air_matched <- subset(data_c_air, vessel_number %in% ves)
data_c_co2_matched <- subset(data_c_co2, vessel_number %in% ves)
data_c_matched <- rbind(data_c_air_matched,data_c_co2_matched)
data.frame(table(data_c_matched$state))   # number of air and co2 should be the same
data_matched <- rbind(data_a,data_v,data_c_matched)

data_matched$id <- factor(data_matched$id)
# data_matched$trt <- factor(data_matched$trt)
data_matched$vessel_number <- factor(data_matched$vessel_number)
data_matched$bolus_number <- factor(data_matched$bolus_number)
data_matched$state <- factor(data_matched$state,levels=c("air","co2"))
data_matched$depth <- factor(data_matched$depth)
nrow(data)

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegressionPCA <- function (df,x_left,x_right) {
  require(ggplot2)
  ggplot(df,aes_string(names(df)[1],names(df)[2]))+geom_point()+
    geom_line(aes(x=air,y=cu),colour="darkgrey",alpha="0.5")+
    geom_line(aes(x=air,y=cl),colour="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cl,ymax=cu),fill="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cu,ymax=cl),fill="darkgrey",alpha="0.5")+
    stat_smooth(data=subset(df,air>x_left&air<x_right),
                method = f, se = TRUE, colour = 'red', size=1.5, formula=y~x)
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
# rat #74 excluded due to imaging depth shift
data <- subset(data,id!=74)
data_matched <- subset(data_matched,id!=74)
airdata_wt_ntx<-subset(data,data$state=="air"&data$gntrt=="wt&ntx")
airdata_tg_ntx<-subset(data,data$state=="air"&data$gntrt=="tg&ntx")
cdata_wt_ntx<-subset(data,data$state=="co2"&data$gntrt=="wt&ntx")
cdata_tg_ntx<-subset(data,data$state=="co2"&data$gntrt=="tg&ntx")
airdata_wt_tx<-subset(data,data$state=="air"&data$gntrt=="wt&tx")
airdata_tg_tx<-subset(data,data$state=="air"&data$gntrt=="tg&tx")
cdata_wt_tx<-subset(data,data$state=="co2"&data$gntrt=="wt&tx")
cdata_tg_tx<-subset(data,data$state=="co2"&data$gntrt=="tg&tx")

# count of each rats and vessels in each group
length(unique(airdata_wt_ntx$id))
length(unique(airdata_wt_ntx[airdata_wt_ntx$type=="a",]$vessel_number))
length(unique(airdata_wt_ntx[airdata_wt_ntx$type=="c",]$vessel_number))
length(unique(airdata_wt_ntx[airdata_wt_ntx$type=="v",]$vessel_number))
length(unique(airdata_tg_ntx$id))
length(unique(airdata_tg_ntx[airdata_tg_ntx$type=="a",]$vessel_number))
length(unique(airdata_tg_ntx[airdata_tg_ntx$type=="c",]$vessel_number))
length(unique(airdata_tg_ntx[airdata_tg_ntx$type=="v",]$vessel_number))
length(unique(airdata_wt_tx$id))
length(unique(airdata_wt_tx[airdata_wt_tx$type=="a",]$vessel_number))
length(unique(airdata_wt_tx[airdata_wt_tx$type=="c",]$vessel_number))
length(unique(airdata_wt_tx[airdata_wt_tx$type=="v",]$vessel_number))
length(unique(airdata_tg_tx$id))
length(unique(airdata_tg_tx[airdata_tg_tx$type=="a",]$vessel_number))
length(unique(airdata_tg_tx[airdata_tg_tx$type=="c",]$vessel_number))
length(unique(airdata_tg_tx[airdata_tg_tx$type=="v",]$vessel_number))

qplot(airdata_wt_ntx$vessel_number,airdata_wt_ntx$pt,colour = airdata_wt_ntx$id)
qplot(airdata_tg_ntx$vessel_number,airdata_tg_ntx$pt,colour = airdata_tg_ntx$id)
qplot(cdata_wt_ntx$vessel_number,cdata_wt_ntx$pt,colour = cdata_wt_ntx$id)
qplot(cdata_tg_ntx$vessel_number,cdata_tg_ntx$pt,colour = cdata_tg_ntx$id)
data_summary_air_wt_ntx <- summarySE(airdata_wt_ntx, measurevar="pt", groupvars=c("depth"))
data_summary_air_wt_ntx
data_summary_air_tg_ntx <- summarySE(airdata_tg_ntx, measurevar="pt", groupvars=c("depth"))
data_summary_air_tg_ntx
data_summary_c_wt_ntx <- summarySE(cdata_wt_ntx, measurevar="pt", groupvars=c("depth"))
data_summary_c_wt_ntx
data_summary_c_tg_ntx <- summarySE(cdata_tg_ntx, measurevar="pt", groupvars=c("depth"))
data_summary_c_tg_ntx

## jitter plots of air and co2 with lines connecting each vessel
data_matched_ntx <- subset(data_matched,trt=="ntx")
data_matched_ntx_a <- subset(data_matched_ntx,state=="air")
data_matched_ntx_c <- subset(data_matched_ntx,state=="co2")
# row index of the measurements with repeating vessel_number
rep_ves_num <- data_matched_ntx_a[duplicated(data_matched_ntx_a$vessel_number),]$vessel_number
ind_air <- which(duplicated(data_matched_ntx_a$vessel_number))
ind_co2 <- which(duplicated(data_matched_ntx_c$vessel_number))
# new set of vessel numbers for repeating vessels
temp_ves_num <- c(seq(max(as.numeric(paste(data_matched_ntx$vessel_number))),
                      max(as.numeric(paste(data_matched_ntx$vessel_number)))-1+
                      length(data_matched_ntx[duplicated(data_matched_ntx$vessel_number)&data_matched_ntx$state=="air",]$vessel_number)))
# assigning the temp_ves_num to the repeating rows
data_matched_ntx_a$vessel_number <- as.numeric(paste(data_matched_ntx_a$vessel_number))
data_matched_ntx_c$vessel_number <- as.numeric(paste(data_matched_ntx_c$vessel_number))
data_matched_ntx_a[ind_air,'vessel_number'] <- temp_ves_num
data_matched_ntx_c[ind_co2,'vessel_number'] <- temp_ves_num
data_matched_ntx <- rbind(data_matched_ntx_a,data_matched_ntx_c)
data_matched_ntx <- data_matched_ntx[order(data_matched_ntx$id,
                                           data_matched_ntx$vessel_number),]
data_matched_ntx_wt <- subset(data_matched_ntx,gn=="wt")
data_matched_ntx_tg <- subset(data_matched_ntx,gn=="tg")
data_matched_ntx_wt_slp <- subset(data_matched_ntx,gn=="wt"&slope!=0&type=="a")
data_matched_ntx_tg_slp <- subset(data_matched_ntx,gn=="tg"&slope!=0&type=="a")
wt_ALL <- ggplot(data_matched_ntx_wt,aes(x=state,y=pt,group=vessel_number
                                         ,color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                    labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-nTg-ALL")
wt_ALL               
tg_ALL <- ggplot(data_matched_ntx_tg,aes(x=state,y=pt,group=vessel_number,
                                     color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                     labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-Tg-ALL")
tg_ALL               
wt_slp_ALL <- ggplot(data_matched_ntx_wt_slp,aes(x=state,y=slope,group=vessel_number
                                         ,color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                     labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-nTg-slope-ALL")+ylim(0,1.3)
wt_slp_ALL               
tg_slp_ALL <- ggplot(data_matched_ntx_tg_slp,aes(x=state,y=slope,group=vessel_number,
                                         color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                     labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-Tg-slope-ALL")+ylim(0,1.3)
tg_slp_ALL
wt_slp_ALL <- ggplot(data_matched_ntx_wt_slp,aes(x=state,y=AUC,group=vessel_number
                                                 ,color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                     labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-nTg-AUC-ALL")+ylim(0,6000000)
wt_slp_ALL               
tg_slp_ALL <- ggplot(data_matched_ntx_tg_slp,aes(x=state,y=AUC,group=vessel_number,
                                                 color=type))+
  scale_color_manual(name="Vessel Types",values=c("red","green","blue"),
                     labels=c("Arterioles","Capillaries","Venules"))+
  geom_line()+ggtitle("nTx-Tg-AUC-ALL")+ylim(0,6000000)
tg_slp_ALL

data_matched_ntx_wt_slp <- subset(data_matched_ntx,gn=="wt"&slope!=0)
data_matched_ntx_tg_slp <- subset(data_matched_ntx,gn=="tg"&slope!=0)
data_matched_ntx_slp <- rbind(data_matched_ntx_tg_slp,data_matched_ntx_wt_slp)

summary(lme(pt ~ slope*state, random = list(~1|id), na.action = na.omit, data_matched_ntx_slp))
summary(lme(pt ~ slope, random = list(~1|id, ~1|state), na.action = na.omit, subset(data_matched_ntx_slp,type=="c")))
summary(lme(slope ~ state, random = list(~1|id), na.action = na.omit, data=data_matched_ntx_slp))
# NO EFFECT         Value  Std.Error  DF   t-value p-value
# stateco2    -0.02679036 0.02964282  60 -0.903772  0.3697 !! No change on slope from air to co2 !!

summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,type=="a")))
summary(lme(AUC ~ gn*state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,type=="a")))
# NO EFFECT         Value  Std.Error  DF   t-value p-value
# stateco2    -0.02679036 0.02964282  60 -0.903772  0.3697 !! No change on slope from air to co2 !!

summary(lme(AUC ~ gn*state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,type=="a")))



airdata <- subset(data_matched_ntx,state=="air"&slope!=0)
cdata <- subset(data_matched_ntx,state=="co2"&slope!=0)
delta_slope = (cdata$slope - airdata$slope)*100/airdata$slope
delta_AUC = (cdata$AUC - airdata$AUC)*100/airdata$AUC
# airdata with slope of co2 and delta_slope concatenated
airdata_co2 <- airdata
colnames(airdata_co2)[10:13] <- c("AUC_air","slope_cl_air","slope_air","slope_cu_air")
target_AUC <- which(names(airdata_co2) == "AUC_air")
target_slope <- which(names(airdata_co2) == "slope_cu_air")
airdata_co2 <- cbind(airdata_co2[,1:target_AUC,drop=F],cdata[,10],delta_AUC,
                     airdata_co2[,(target_AUC+1):target_slope,drop=F],cdata[,11:13],delta_slope,
                     airdata_co2[,(target_slope+1):length(airdata_co2),drop=F])
colnames(airdata_co2)[11] <- c("AUC_co2")
colnames(airdata_co2)[16:19] <- c("slope_cl_co2","slope_co2","slope_cu_co2","delta_slope")
airdata_co2 <- airdata_co2[order(rev(airdata_co2$gn),airdata_co2$id,airdata_co2$vessel_number),]
airdata_co2$gn <- factor(airdata_co2$gn,c("wt","tg"))
airdata_co2_all <- airdata_co2
tmp <- as.character(airdata_co2_all$type)
tmp <- c("All")
airdata_co2_all$type <- factor(tmp)
airdata_co2_all <- rbind(airdata_co2,airdata_co2_all)
airdata_co2_all$type = with(airdata_co2_all, factor(type, levels = rev(levels(type))))
# errors in data_summary_delta not correct since propagation of error not done
data_summary_delta <- summarySE(airdata_co2_all, measurevar="delta_slope", groupvars=c("type","gn"))
data_summary_air <- summarySE(airdata_co2_all, measurevar="slope_air", groupvars=c("type","gn"))
data_summary_co2 <- summarySE(airdata_co2_all, measurevar="slope_co2", groupvars=c("type","gn"))
data_summary_delta_AUC <- summarySE(airdata_co2_all, measurevar="delta_AUC", groupvars=c("type","gn"))
data_summary_air_AUC <- summarySE(airdata_co2_all, measurevar="AUC_air", groupvars=c("type","gn"))
data_summary_co2_AUC <- summarySE(airdata_co2_all, measurevar="AUC_co2", groupvars=c("type","gn"))
# calculating correct erros for data_summary_delta following error propagation rules
# Doing things manually because I can
# covariance for each group
data_summary_delta_corr <- data_summary_delta
cov_summary <- c(cov(subset(airdata_co2_all,type=="All"&gn=="wt")$slope_air,subset(airdata_co2_all,type=="All"&gn=="wt")$slope_co2),
                 cov(subset(airdata_co2_all,type=="All"&gn=="tg")$slope_air,subset(airdata_co2_all,type=="All"&gn=="tg")$slope_co2),
                 cov(subset(airdata_co2_all,type=="v"&gn=="wt")$slope_air,subset(airdata_co2_all,type=="v"&gn=="wt")$slope_co2),
                 cov(subset(airdata_co2_all,type=="v"&gn=="tg")$slope_air,subset(airdata_co2_all,type=="v"&gn=="tg")$slope_co2),
                 cov(subset(airdata_co2_all,type=="c"&gn=="wt")$slope_air,subset(airdata_co2_all,type=="c"&gn=="wt")$slope_co2),
                 cov(subset(airdata_co2_all,type=="c"&gn=="tg")$slope_air,subset(airdata_co2_all,type=="c"&gn=="tg")$slope_co2),
                 cov(subset(airdata_co2_all,type=="a"&gn=="wt")$slope_air,subset(airdata_co2_all,type=="a"&gn=="wt")$slope_co2),
                 cov(subset(airdata_co2_all,type=="a"&gn=="tg")$slope_air,subset(airdata_co2_all,type=="a"&gn=="tg")$slope_co2))
data_summary_delta_AUC_corr <- data_summary_delta_AUC
cov_summary_AUC <- c(cov(subset(airdata_co2_all,type=="All"&gn=="wt")$AUC_air,subset(airdata_co2_all,type=="All"&gn=="wt")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="All"&gn=="tg")$AUC_air,subset(airdata_co2_all,type=="All"&gn=="tg")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="v"&gn=="wt")$AUC_air,subset(airdata_co2_all,type=="v"&gn=="wt")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="v"&gn=="tg")$AUC_air,subset(airdata_co2_all,type=="v"&gn=="tg")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="c"&gn=="wt")$AUC_air,subset(airdata_co2_all,type=="c"&gn=="wt")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="c"&gn=="tg")$AUC_air,subset(airdata_co2_all,type=="c"&gn=="tg")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="a"&gn=="wt")$AUC_air,subset(airdata_co2_all,type=="a"&gn=="wt")$AUC_co2),
                 cov(subset(airdata_co2_all,type=="a"&gn=="tg")$AUC_air,subset(airdata_co2_all,type=="a"&gn=="tg")$AUC_co2))
# standard deviation following error propagation
data_summary_delta_corr$sd[1] <- sqrt((data_summary_co2$sd[1]/data_summary_co2$slope_co2[1])^2
                                +(data_summary_air$sd[1]/data_summary_air$slope_air[1])^2
                                -2*cov_summary[1]/(data_summary_co2$slope_co2[1]*data_summary_air$slope_air[1]))*
                                abs(data_summary_delta_corr$delta_slope[1])
data_summary_delta_corr$sd[2] <- sqrt((data_summary_co2$sd[2]/data_summary_co2$slope_co2[2])^2
                                +(data_summary_air$sd[2]/data_summary_air$slope_air[2])^2
                                -2*cov_summary[2]/(data_summary_co2$slope_co2[2]*data_summary_air$slope_air[2]))*
                                abs(data_summary_delta_corr$delta_slope[2])
data_summary_delta_corr$sd[3] <- sqrt((data_summary_co2$sd[3]/data_summary_co2$slope_co2[3])^2
                                +(data_summary_air$sd[3]/data_summary_air$slope_air[3])^2
                                -2*cov_summary[3]/(data_summary_co2$slope_co2[3]*data_summary_air$slope_air[3]))*
                                abs(data_summary_delta_corr$delta_slope[3])
data_summary_delta_corr$sd[4] <- sqrt((data_summary_co2$sd[4]/data_summary_co2$slope_co2[4])^2
                                +(data_summary_air$sd[4]/data_summary_air$slope_air[4])^2
                                -2*cov_summary[4]/(data_summary_co2$slope_co2[4]*data_summary_air$slope_air[4]))*
                                abs(data_summary_delta_corr$delta_slope[4])
data_summary_delta_corr$sd[5] <- sqrt((data_summary_co2$sd[5]/data_summary_co2$slope_co2[5])^2
                                +(data_summary_air$sd[5]/data_summary_air$slope_air[5])^2
                                -2*cov_summary[5]/(data_summary_co2$slope_co2[5]*data_summary_air$slope_air[5]))*
                                abs(data_summary_delta_corr$delta_slope[5])
data_summary_delta_corr$sd[6] <- sqrt((data_summary_co2$sd[6]/data_summary_co2$slope_co2[6])^2
                                +(data_summary_air$sd[6]/data_summary_air$slope_air[6])^2
                                -2*cov_summary[6]/(data_summary_co2$slope_co2[6]*data_summary_air$slope_air[6]))*
                                abs(data_summary_delta_corr$delta_slope[6])
data_summary_delta_corr$sd[7] <- sqrt((data_summary_co2$sd[7]/data_summary_co2$slope_co2[7])^2
                                +(data_summary_air$sd[7]/data_summary_air$slope_air[7])^2
                                -2*cov_summary[7]/(data_summary_co2$slope_co2[7]*data_summary_air$slope_air[7]))*
                                abs(data_summary_delta_corr$delta_slope[7])
data_summary_delta_corr$sd[8] <- sqrt((data_summary_co2$sd[8]/data_summary_co2$slope_co2[8])^2
                                +(data_summary_air$sd[8]/data_summary_air$slope_air[8])^2
                                -2*cov_summary[8]/(data_summary_co2$slope_co2[8]*data_summary_air$slope_air[8]))*
                                abs(data_summary_delta_corr$delta_slope[8])
data_summary_delta_corr$se <- data_summary_delta_corr$sd / sqrt(data_summary_delta_corr$N)  # Calculate new standard error of the mean
ciMult <- qt(0.95/2 + .5, data_summary_delta_corr$N-1)
data_summary_delta_corr$ci <- data_summary_delta_corr$se * ciMult

data_summary_delta_AUC_corr$sd[1] <- sqrt((data_summary_co2_AUC$sd[1]/data_summary_co2_AUC$AUC_co2[1])^2
                                      +(data_summary_air_AUC$sd[1]/data_summary_air_AUC$AUC_air[1])^2
                                      -2*cov_summary_AUC[1]/(data_summary_co2_AUC$AUC_co2[1]*data_summary_air_AUC$AUC_air[1]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[1])
data_summary_delta_AUC_corr$sd[2] <- sqrt((data_summary_co2_AUC$sd[2]/data_summary_co2_AUC$AUC_co2[2])^2
                                      +(data_summary_air_AUC$sd[2]/data_summary_air_AUC$AUC_air[2])^2
                                      -2*cov_summary_AUC[2]/(data_summary_co2_AUC$AUC_co2[2]*data_summary_air_AUC$AUC_air[2]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[2])
data_summary_delta_AUC_corr$sd[3] <- sqrt((data_summary_co2_AUC$sd[3]/data_summary_co2_AUC$AUC_co2[3])^2
                                      +(data_summary_air_AUC$sd[3]/data_summary_air_AUC$AUC_air[3])^2
                                      -2*cov_summary_AUC[3]/(data_summary_co2_AUC$AUC_co2[3]*data_summary_air_AUC$AUC_air[3]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[3])
data_summary_delta_AUC_corr$sd[4] <- sqrt((data_summary_co2_AUC$sd[4]/data_summary_co2_AUC$AUC_co2[4])^2
                                      +(data_summary_air_AUC$sd[4]/data_summary_air_AUC$AUC_air[4])^2
                                      -2*cov_summary_AUC[4]/(data_summary_co2_AUC$AUC_co2[4]*data_summary_air_AUC$AUC_air[4]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[4])
data_summary_delta_AUC_corr$sd[5] <- sqrt((data_summary_co2_AUC$sd[5]/data_summary_co2_AUC$AUC_co2[5])^2
                                      +(data_summary_air_AUC$sd[5]/data_summary_air_AUC$AUC_air[5])^2
                                      -2*cov_summary_AUC[5]/(data_summary_co2_AUC$AUC_co2[5]*data_summary_air_AUC$AUC_air[5]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[5])
data_summary_delta_AUC_corr$sd[6] <- sqrt((data_summary_co2_AUC$sd[6]/data_summary_co2_AUC$AUC_co2[6])^2
                                      +(data_summary_air_AUC$sd[6]/data_summary_air_AUC$AUC_air[6])^2
                                      -2*cov_summary_AUC[6]/(data_summary_co2_AUC$AUC_co2[6]*data_summary_air_AUC$AUC_air[6]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[6])
data_summary_delta_AUC_corr$sd[7] <- sqrt((data_summary_co2_AUC$sd[7]/data_summary_co2_AUC$AUC_co2[7])^2
                                      +(data_summary_air_AUC$sd[7]/data_summary_air_AUC$AUC_air[7])^2
                                      -2*cov_summary_AUC[7]/(data_summary_co2_AUC$AUC_co2[7]*data_summary_air_AUC$AUC_air[7]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[7])
data_summary_delta_AUC_corr$sd[8] <- sqrt((data_summary_co2_AUC$sd[8]/data_summary_co2_AUC$AUC_co2[8])^2
                                      +(data_summary_air_AUC$sd[8]/data_summary_air_AUC$AUC_air[8])^2
                                      -2*cov_summary_AUC[8]/(data_summary_co2_AUC$AUC_co2[8]*data_summary_air_AUC$AUC_air[8]))*
                                      abs(data_summary_delta_AUC_corr$delta_AUC[8])
data_summary_delta_AUC_corr$se <- data_summary_delta_AUC_corr$sd / sqrt(data_summary_delta_AUC_corr$N)  # Calculate new standard error of the mean
ciMult <- qt(0.95/2 + .5, data_summary_delta_AUC_corr$N-1)
data_summary_delta_AUC_corr$ci <- data_summary_delta_AUC_corr$se * ciMult






summary(lme(delta_slope ~ gn, random = ~1 | id, na.action = na.omit, subset(airdata_co2_all,type=="a")))
#                   Value Std.Error  DF t-value  p-value
#    gntg        63.29230  51.16826  12 1.236944  0.2398
summary(lme(delta_AUC ~ gn, random = ~1 | id/bolus_number, na.action = na.omit, subset(airdata_co2_all,type=="a")))
#                   Value Std.Error  DF t-value  p-value
#    gntg         61.28553  27.40271  12  2.236477  0.0451
airdata_co2_all_all <- subset(airdata_co2_all,type=="All")
wilcox.test(delta_slope~gn,data=airdata_co2_all_all)
airdata_co2_all_a <- subset(airdata_co2_all,type=="a")
wilcox.test(delta_slope~gn,data=airdata_co2_all_a)
airdata_co2_all_c <- subset(airdata_co2_all,type=="c")
wilcox.test(delta_slope~gn,data=airdata_co2_all_c)
airdata_co2_all_v <- subset(airdata_co2_all,type=="v")
wilcox.test(delta_slope~gn,data=airdata_co2_all_v)
airdata_co2_all_all <- subset(airdata_co2_all,type=="All")
wilcox.test(delta_AUC~gn,data=airdata_co2_all_all)
airdata_co2_all_a <- subset(airdata_co2_all,type=="a")
wilcox.test(delta_AUC~gn,data=airdata_co2_all_a)
airdata_co2_all_c <- subset(airdata_co2_all,type=="c")
wilcox.test(delta_AUC~gn,data=airdata_co2_all_c)
airdata_co2_all_v <- subset(airdata_co2_all,type=="v")
wilcox.test(delta_AUC~gn,data=airdata_co2_all_v)

# Confidence interval multiplier for standard error
# Calculate t-statistic for confidence interval: 
# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1

aa<-ggplot(data_summ, aes(x=type, y=delta_slope)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=gn,ymin=delta_slope, ymax=delta_slope+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab(expression(paste(Delta, "TTP (%)"))) + ylim(-15,120)+
  scale_fill_manual(name="Genotype",values=c("green","red"),labels=c("nTg","Tg"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))
aa

airdata_ntx<-subset(data_matched,data_matched$state=="air"&data_matched$trt=="ntx")
cdata_ntx<-subset(data_matched,data_matched$state=="co2"&data_matched$trt=="ntx")
delta_ntx_pt <- airdata_ntx$pt - cdata_ntx$pt 
airdata_ntx$delta_pt <- delta_ntx_pt
airdata_ntx <- airdata_ntx[order(rev(airdata_ntx$gn),airdata_ntx$id,airdata_ntx$vessel_number),] #
airdata_ntx <- subset(airdata_ntx,id!=74)
airdata_ntx$gn <- factor(airdata_ntx$gn,c("wt","tg"))
airdata_ntx_all <- airdata_ntx
tmp <- as.character(airdata_ntx_all$type)
tmp <- "All"
airdata_ntx_all$type <- factor(tmp)
airdata_ntx_all <- rbind(airdata_ntx,airdata_ntx_all)
airdata_ntx_all$type = with(airdata_ntx_all, factor(type, levels = rev(levels(type))))
data_summary <- summarySE(airdata_ntx_all, measurevar="delta_pt", groupvars=c("type","gn"))
data_summ <- data_summary
data_summ
# TTP bar
aa<-ggplot(data_summ, aes(x=type, y=delta_pt)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=gn,ymin=delta_pt, ymax=delta_pt+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab("??TTP (s)") + ylim(0,2.5)+
  scale_fill_manual(name="Genotype",values=c("green","red"),labels=c("nTg","Tg"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))






# 
# by_vessel_n <-split(data$pt,data$vessel_number)
# 
# airdata$TP_rest_centered <- NA
# cdata$TP_CO2_centered <- NA
# 
# qplot(data$vessel_number,data$pt,color=data$state)
# delta_wt_ntx = cdata_wt_ntx$pt - airdata_wt_ntx$pt
# delta_tg_ntx = cdata_wt_ntx$pt - airdata_wt_ntx$pt
# 
# airdata_wt_ntx$delta_pt <- delta_wt_ntx
# airdata_tg_ntx$delta_pt <- delta_tg_ntx
# vessel_number = 1:97
# aa<-qplot(airdata$vessel_number,airdata$delta_pt,grid,colour = airdata$id)
# aa+labs(y="delta_pt [s]")+labs(x="vessel_number")+labs(title="TTPco2 - TTPair")

airdata <- subset(data_matched,state=="air")
cdata <- subset(data_matched,state=="co2")

delta_pt = cdata$pt - airdata$pt
airdata$delta_pt <- delta_pt

bb<-ggplot(data = airdata,aes(x=gntrt,y=delta_pt,colour=gntrt))+geom_boxplot()
bb+labs(y=expression(paste(Delta,"TTP (s)")))+
  labs(x="Genotype & Treatment")+
  theme(axis.text=element_text(size=13),
        text=element_text(size=18),
        legend.text=element_text(size=15))

airdata_a <- subset(airdata,type=="a")
cdata_a <- subset(cdata,type=="a")
delta_a_pt <- cdata_a$pt - airdata_a$pt
airdata_a$delta_pt <- delta_a_pt
cc<-ggplot(data = airdata_a,aes(x=gntrt,y=delta_pt,colour=gntrt))+geom_boxplot()
cc+labs(y=expression(paste(Delta,"TTP (s)")))+
  labs(x="Genotype & Treatment")+
  theme(axis.text=element_text(size=13),
        text=element_text(size=18),
        legend.text=element_text(size=15))
## point plots to see if separation of datapoints by vessel type is visibly clear (it's not)
# airdata_ntx_wt_all <- airdata_ntx_wt
# tmp <- as.character(airdata_ntx_wt_all$type)
# tmp <- "All"
# airdata_ntx_wt_all$type <- factor(tmp)
# airdata_ntx_wt_all <- rbind(airdata_ntx_wt,airdata_ntx_wt_all)
# airdata_ntx_wt_all$type = with(airdata_ntx_wt_all, factor(type, levels = rev(levels(type))))
# airdata_ntx_wt_all$id = factor(airdata_ntx_wt_all$id)
# bb<-ggplot(data = airdata_ntx_wt_all,aes(x=type,y=pt))+geom_point()
# bb+labs(y=expression(paste("TTP (s)")))+
#   labs(x="Vessel Type")+
#   theme(axis.text=element_text(size=15),
#         text=element_text(size=20),
#         legend.text=element_text(size=18))+
#   ggtitle(paste("Air Breathing - nTg"))

## jitter plot by state (air & co2)
data_wt_ntx<-subset(data,data$gntrt=="wt&ntx")
data_wt_ntx_art<-subset(data_wt_ntx,type=="a")
data_wt_ntx_cap<-subset(data_wt_ntx,type=="c")
data_wt_ntx_ven<-subset(data_wt_ntx,type=="v")
data_tg_ntx<-subset(data,data$gntrt=="tg&ntx")
data_tg_ntx_art<-subset(data_tg_ntx,type=="a")
data_tg_ntx_cap<-subset(data_tg_ntx,type=="c")
data_tg_ntx_ven<-subset(data_tg_ntx,type=="v")
vessel_names <- list(
  'a'="Arterioles",
  'c'="Capillaries",
  'v'="Venules")
vessel_labeller <- function(variable,value){
  return(vessel_names[value])
}

data_wt_ntx_summ <- summarySE(data_wt_ntx, measurevar="pt", groupvars=c("state","type"))
data_tg_ntx_summ <- summarySE(data_tg_ntx, measurevar="pt", groupvars=c("state","type"))
data_wt_ntx_art_summ <- summarySE(data_wt_ntx_art, measurevar="pt", groupvars=c("state"))
jitter_wt_ntx <- ggplot(data_wt_ntx,aes(x=state,y=pt,color=state))+ylim(5,12)+labs(x="",y="TTP (s)")+
  geom_jitter(position=position_jitter(w=0.1,h=0.1))+
  layer(data=data_wt_ntx_summ,mapping=aes(x=state,y=pt),geom="point",size=5,color="red")+
  facet_grid(. ~ type,scales="free",space="free",labeller=vessel_labeller)+
  scale_x_discrete(labels=c("Air","CO2"))+
  ggtitle(paste("nTg"))+
  theme(axis.text=element_text(size=16,colour="black"),
        text=element_text(size=24),
        legend.text=element_text(size=19))
jitter_wt_ntx

data_tg_ntx_art_summ <- summarySE(data_tg_ntx_art, measurevar="pt", groupvars=c("state"))
jitter_tg_ntx <- ggplot(data_tg_ntx,aes(x=state,y=pt,color=state))+ylim(5,12)+labs(x="",y="TTP (s)")+
  geom_jitter(position=position_jitter(w=0.1,h=0.1))+
  layer(data=data_tg_ntx_summ,mapping=aes(x=state,y=pt),geom="point",size=5,color="red")+
  facet_grid(. ~ type,scales="free",space="free",labeller=vessel_labeller)+
  scale_x_discrete(labels=c("Air","CO2"))+
  ggtitle(paste("Tg"))+
  theme(axis.text=element_text(size=16,colour="black"),
        text=element_text(size=24),
        legend.text=element_text(size=19))
jitter_tg_ntx
# 
# p + geom_jitter(position = position_jitter(w = 0.1, h = 0.1)) + 
#   layer(data = a,mapping = aes(x = female, y = y.mean), geom = "point", size = 5, color = "red") +
#   facet_grid(. ~ subject, scales = "free", space = "free") + theme_bw()

## density plots by vessel types
data_ntx = subset(data_matched,trt=="ntx")
airdata_ntx_wt <- subset(data_ntx,state=="air"&gn=="wt")
cdata_ntx_wt <- subset(data_ntx,state=="co2"&gn=="wt")
airdata_ntx_tg <- subset(data_ntx,state=="air"&gn=="tg")
cdata_ntx_tg <- subset(data_ntx,state=="co2"&gn=="tg")
data_ntx_wt_comb <- rbind(airdata_ntx_wt,cdata_ntx_wt)
aa<-ggplot(data = data_ntx_wt_comb,aes(x=type,y=pt))+geom_violin(aes(fill=state))
aa<-aa+labs(y=expression(paste("TTP (s)")))+
  labs(x="Vessel Type")+ylim(4,12)+guides(colour=F)+
  theme(axis.text=element_text(size=15,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))+ggtitle(paste("nTg"))+
  scale_x_discrete(labels=c("Arterioles","Capillaries","Venules"))
aa
data_ntx_tg_comb <- rbind(airdata_ntx_tg,cdata_ntx_tg)
bb<-ggplot(data = data_ntx_tg_comb,aes(x=type,y=pt))+geom_violin(aes(fill=state))
bb<-bb+labs(y=expression(paste("TTP (s)")))+
  labs(x="Vessel Type")+ylim(4,12)+guides(colour=F)+
  theme(axis.text=element_text(size=15,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))+ggtitle(paste("Tg"))+
  scale_x_discrete(labels=c("Arterioles","Capillaries","Venules"))
bb
grid.arrange(aa,bb,ncol=1,nrow=2)
dev.copy(jpeg,"TTP_homogenization_nTgTg.jpeg")
dev.off()
sd_wt_ntx_air_a <- sd(airdata_wt_ntx[airdata_wt_ntx$type=="a",]$pt)
sd_wt_ntx_co2_a <- sd(cdata_wt_ntx[cdata_wt_ntx$type=="a",]$pt)
sd_wt_ntx_air_c <- sd(airdata_wt_ntx[airdata_wt_ntx$type=="c",]$pt)
sd_wt_ntx_co2_c <- sd(cdata_wt_ntx[cdata_wt_ntx$type=="c",]$pt)
sd_wt_ntx_air_v <- sd(airdata_wt_ntx[airdata_wt_ntx$type=="v",]$pt)
sd_wt_ntx_co2_v <- sd(cdata_wt_ntx[cdata_wt_ntx$type=="v",]$pt)
sd_tg_ntx_air_a <- sd(airdata_tg_ntx[airdata_tg_ntx$type=="a",]$pt)
sd_tg_ntx_co2_a <- sd(cdata_tg_ntx[cdata_tg_ntx$type=="a",]$pt)
sd_tg_ntx_air_c <- sd(airdata_tg_ntx[airdata_tg_ntx$type=="c",]$pt)
sd_tg_ntx_co2_c <- sd(cdata_tg_ntx[cdata_tg_ntx$type=="c",]$pt)
sd_tg_ntx_air_v <- sd(airdata_tg_ntx[airdata_tg_ntx$type=="v",]$pt)
sd_tg_ntx_co2_v <- sd(cdata_tg_ntx[cdata_tg_ntx$type=="v",]$pt)
mean_wt_ntx_air_a <- mean(airdata_wt_ntx[airdata_wt_ntx$type=="a",]$pt)
mean_wt_ntx_co2_a <- mean(cdata_wt_ntx[cdata_wt_ntx$type=="a",]$pt)
mean_wt_ntx_air_c <- mean(airdata_wt_ntx[airdata_wt_ntx$type=="c",]$pt)
mean_wt_ntx_co2_c <- mean(cdata_wt_ntx[cdata_wt_ntx$type=="c",]$pt)
mean_wt_ntx_air_v <- mean(airdata_wt_ntx[airdata_wt_ntx$type=="v",]$pt)
mean_wt_ntx_co2_v <- mean(cdata_wt_ntx[cdata_wt_ntx$type=="v",]$pt)
mean_tg_ntx_air_a <- mean(airdata_tg_ntx[airdata_tg_ntx$type=="a",]$pt)
mean_tg_ntx_co2_a <- mean(cdata_tg_ntx[cdata_tg_ntx$type=="a",]$pt)
mean_tg_ntx_air_c <- mean(airdata_tg_ntx[airdata_tg_ntx$type=="c",]$pt)
mean_tg_ntx_co2_c <- mean(cdata_tg_ntx[cdata_tg_ntx$type=="c",]$pt)
mean_tg_ntx_air_v <- mean(airdata_tg_ntx[airdata_tg_ntx$type=="v",]$pt)
mean_tg_ntx_co2_v <- mean(cdata_tg_ntx[cdata_tg_ntx$type=="v",]$pt)
cov_wt_ntx_air_a <- sd_wt_ntx_air_a/mean_wt_ntx_air_a
cov_wt_ntx_co2_a <- sd_wt_ntx_co2_a/mean_wt_ntx_co2_a
cov_wt_ntx_air_c <- sd_wt_ntx_air_c/mean_wt_ntx_air_c
cov_wt_ntx_co2_c <- sd_wt_ntx_co2_c/mean_wt_ntx_co2_c
cov_wt_ntx_air_v <- sd_wt_ntx_air_v/mean_wt_ntx_air_v
cov_wt_ntx_co2_v <- sd_wt_ntx_co2_v/mean_wt_ntx_co2_v
cov_tg_ntx_air_a <- sd_tg_ntx_air_a/mean_tg_ntx_air_a
cov_tg_ntx_co2_a <- sd_tg_ntx_co2_a/mean_tg_ntx_co2_a
cov_tg_ntx_air_c <- sd_tg_ntx_air_c/mean_tg_ntx_air_c
cov_tg_ntx_co2_c <- sd_tg_ntx_co2_c/mean_tg_ntx_co2_c
cov_tg_ntx_air_v <- sd_tg_ntx_air_v/mean_tg_ntx_air_v
cov_tg_ntx_co2_v <- sd_tg_ntx_co2_v/mean_tg_ntx_co2_v
sd_wt <- c(sd_wt_ntx_air_a,sd_wt_ntx_co2_a,sd_wt_ntx_air_c,sd_wt_ntx_co2_c,
           sd_wt_ntx_air_v,sd_wt_ntx_co2_v)
sd_tg <- c(sd_tg_ntx_air_a,sd_tg_ntx_co2_a,sd_tg_ntx_air_c,sd_tg_ntx_co2_c,
           sd_tg_ntx_air_v,sd_tg_ntx_co2_v)
sd_all <- t(data.frame(sd_wt,sd_tg))
cov_wt <- c(cov_wt_ntx_air_a,cov_wt_ntx_co2_a,cov_wt_ntx_air_c,cov_wt_ntx_co2_c,
           cov_wt_ntx_air_v,cov_wt_ntx_co2_v)
cov_tg <- c(cov_tg_ntx_air_a,cov_tg_ntx_co2_a,cov_tg_ntx_air_c,cov_tg_ntx_co2_c,
           cov_tg_ntx_air_v,cov_tg_ntx_co2_v)
cov_all <- t(data.frame(cov_wt,cov_tg))
var.test(airdata_wt_ntx[airdata_wt_ntx$type=="a",]$pt,cdata_wt_ntx[cdata_wt_ntx$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_wt_ntx[airdata_wt_ntx$type=="c",]$pt,cdata_wt_ntx[cdata_wt_ntx$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_wt_ntx[airdata_wt_ntx$type=="v",]$pt,cdata_wt_ntx[cdata_wt_ntx$type=="v",]$pt,alternative=c("greater"))
var.test(airdata_tg_ntx[airdata_tg_ntx$type=="a",]$pt,cdata_tg_ntx[cdata_tg_ntx$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_tg_ntx[airdata_tg_ntx$type=="c",]$pt,cdata_tg_ntx[cdata_tg_ntx$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_tg_ntx[airdata_tg_ntx$type=="v",]$pt,cdata_tg_ntx[cdata_tg_ntx$type=="v",]$pt,alternative=c("greater"))

# FC_60sec <- c(paste(round(pct50_60sec_flow_change,digits=3)," +/- ",round(pct50_60sec_flow_change_std,digits=3)),
#               paste(round(pct75_60sec_flow_change,digits=3)," +/- ",round(pct75_60sec_flow_change_std,digits=3)),
#               paste(round(pct100_60sec_flow_change,digits=3)," +/- ",round(pct100_60sec_flow_change_std,digits=3)))
# FC_60sec_cap <- c(paste(round(pct50_60sec_cap_flow_change,digits=3)," +/- ",round(pct50_60sec_cap_flow_change_std,digits=3)),
#                   paste(round(pct75_60sec_cap_flow_change,digits=3)," +/- ",round(pct75_60sec_cap_flow_change_std,digits=3)),
#                   paste(round(pct100_60sec_cap_flow_change,digits=3)," +/- ",round(pct100_60sec_cap_flow_change_std,digits=3)))
# FC_60sec_pv <- c(paste(round(pct50_60sec_pv_flow_change,digits=3)," +/- ",round(pct50_60sec_pv_flow_change_std,digits=3)),
#                  paste(round(pct75_60sec_pv_flow_change,digits=3)," +/- ",round(pct75_60sec_pv_flow_change_std,digits=3)),
#                  paste(round(pct100_60sec_pv_flow_change,digits=3)," +/- ",round(pct100_60sec_pv_flow_change_std,digits=3)))
# FC_60sec_df <-data.frame(FC_60sec,FC_60sec_cap,FC_60sec_pv)
# setattr(FC_60sec_df,"row.names",c("5% CO2 - 60sec","7.5% CO2 - 60sec","10% CO2 - 60sec"))
# 
# FC_60sec_df.table <- xtable(FC_60sec_df[1:3,])
# print.xtable(FC_60sec_df.table,floating=FALSE,type="latex",file="aaa.tex")


# barplot change in ttp by genotype for ntx only
airdata_ntx<-subset(data_matched,data_matched$state=="air"&data_matched$trt=="ntx")
cdata_ntx<-subset(data_matched,data_matched$state=="co2"&data_matched$trt=="ntx")
delta_ntx_pt <- airdata_ntx$pt - cdata_ntx$pt 
airdata_ntx$delta_pt <- delta_ntx_pt
airdata_ntx <- airdata_ntx[order(rev(airdata_ntx$gn),airdata_ntx$id,airdata_ntx$vessel_number),] #
airdata_ntx <- subset(airdata_ntx,id!=74)
airdata_ntx$gn <- factor(airdata_ntx$gn,c("wt","tg"))
airdata_ntx_all <- airdata_ntx
tmp <- as.character(airdata_ntx_all$type)
tmp <- "All"
airdata_ntx_all$type <- factor(tmp)
airdata_ntx_all <- rbind(airdata_ntx,airdata_ntx_all)
airdata_ntx_all$type = with(airdata_ntx_all, factor(type, levels = rev(levels(type))))
data_summary <- summarySE(airdata_ntx_all, measurevar="delta_pt", groupvars=c("type","gn"))
data_summ <- data_summary
data_summ
# TTP bar
aa<-ggplot(data_summ, aes(x=type, y=delta_pt)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=gn,ymin=delta_pt, ymax=delta_pt+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab("??TTP (s)") + ylim(0,2.5)+
  scale_fill_manual(name="Genotype",values=c("green","red"),labels=c("nTg","Tg"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
      text=element_text(size=20),
      legend.text=element_text(size=15))
label1.df <- data.frame(type = factor(c("All", "a")),
                       delta_pt = c(1.5,1.6))
label2.df <- data.frame(type = factor(c("v")),
                        delta_pt = c(2.1))
data1 <- data.frame(x = c(0.75,0.75,1.25,1.25), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(1.75,1.75,2.25,2.25), y = c(1.5, 1.6, 1.6, 1.5),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(3.75,3.75,4.25,4.25), y = c(2,2.1,2.1,2),
                    type=NA, delta_pt=NA)
aa<-aa + geom_path(data = data1, aes(x = x, y = y),size=1)+
  geom_path(data = data2, aes(x = x, y = y),size=1)+
  geom_path(data = data3, aes(x = x, y = y),size=1)+
  geom_text(data=label1.df,label="**",size=15)+
  geom_text(data=label2.df,label="*",size=15)
aa
dev.copy(jpeg,"??TTP_all_wtVStg.jpeg")
dev.off()
airdata_ntx_all_all <- subset(airdata_ntx_all,type=="All")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_all)
airdata_ntx_all_art <- subset(airdata_ntx_all,type=="a")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_art)
airdata_ntx_all_cap <- subset(airdata_ntx_all,type=="c")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_cap)
airdata_ntx_all_ven <- subset(airdata_ntx_all,type=="v")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_ven)

# barplot change in ttp by genotype and treatment
airdata_all<-subset(data_matched,data_matched$state=="air")
cdata_all<-subset(data_matched,data_matched$state=="co2")
delta_all_pt <- airdata_all$pt - cdata_all$pt 
airdata_all$delta_pt <- delta_all_pt
airdata_all <- airdata_all[order(rev(airdata_all$gntrt),airdata_all$id,airdata_all$vessel_number),] #
airdata_all <- subset(airdata_all,id!=74)
airdata_all$gntrt <- factor(airdata_all$gntrt,c("wt&ntx","tg&ntx","wt&tx","tg&tx"))
airdata_all_all <- airdata_all
tmp <- as.character(airdata_all_all$type)
tmp <- "All"
airdata_all_all$type <- factor(tmp)
airdata_all_all <- rbind(airdata_all,airdata_all_all)
airdata_all_all$type = with(airdata_all_all, factor(type, levels = rev(levels(type))))
airdata_all_ntx <- subset(airdata_all_all,trt=="ntx")
airdata_all_ntx$gn = with(airdata_all_ntx, factor(gn, levels = rev(levels(gn))))
airdata_all_ntx <- subset(airdata_all_all,trt=="ntx")
airdata_all_ntx$gn = with(airdata_all_ntx, factor(gn, levels = rev(levels(gn))))
airdata_all_wt <- subset(airdata_all_all,gn=="wt")
airdata_all_wt$trt = with(airdata_all_wt, factor(trt, levels = rev(levels(trt))))
data_summary <- summarySE(airdata_all_all, measurevar="delta_pt", groupvars=c("type","gntrt"))
data_summ <- data_summary
data_summ
data_summary <- summarySE(airdata_all_ntx, measurevar="delta_pt", groupvars=c("type","gn"))
data_summ <- data_summary
data_summ
data_summary <- summarySE(airdata_all_wt, measurevar="delta_pt", groupvars=c("type","trt"))
data_summ <- data_summary
data_summ
# TTP bar
# 2 bars in each group for ntx ntg vs tg
aa<-ggplot(data_summ, aes(x=type, y=delta_pt)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=gn,ymin=delta_pt, ymax=delta_pt+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab("??TTP (s)") + ylim(0,2.5)+
  scale_fill_manual(name="Genotype",values=c("green","red"),labels=c("nTg","Tg"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))
aa
airdata_all_ntx_all <- subset(airdata_all_ntx,type=="All")
wilcox.test(delta_pt~gn,data=airdata_all_ntx_all)
airdata_all_ntx_art <- subset(airdata_all_ntx,type=="a")
wilcox.test(delta_pt~gn,data=airdata_all_ntx_art)
airdata_all_ntx_cap <- subset(airdata_all_ntx,type=="c")
wilcox.test(delta_pt~gn,data=airdata_all_ntx_cap)
airdata_all_ntx_ven <- subset(airdata_all_ntx,type=="v")
wilcox.test(delta_pt~gn,data=airdata_all_ntx_ven)
label1.df <- data.frame(type = factor(c("All", "a","v")),
                        delta_pt = c(1.5,1.6,2.1))
label2.df <- data.frame(type = factor(c("c")),
                        delta_pt = c(1.5))
data1 <- data.frame(x = c(0.77,0.77,1.23,1.23), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(1.77,1.77,2.23,2.23), y = c(1.5, 1.6, 1.6, 1.5),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(2.77,2.77,3.23,3.23), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data4 <- data.frame(x = c(3.77,3.77,4.23,4.23), y = c(2,2.1,2.1,2),
                    type=NA, delta_pt=NA)
aa<-aa + geom_path(data = data1, aes(x = x, y = y),size=1)+
  geom_path(data = data2, aes(x = x, y = y),size=1)+
  geom_path(data = data3, aes(x = x, y = y),size=1)+
  geom_path(data = data4, aes(x = x, y = y),size=1)+
  geom_text(data=label1.df,label="**",size=15)+
  geom_text(data=label2.df,label="*",size=15)
aa
# 2 bars in each group for ntg untreated vs treated
aa<-ggplot(data_summ, aes(x=type, y=delta_pt)) +
  geom_bar(aes(fill=trt),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=trt,ymin=delta_pt, ymax=delta_pt+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab("??TTP (s)") + ylim(0,2.5)+
  scale_fill_manual(name="Treatment",values=c("white","black"),labels=c("no LNAME","LNAME"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))
aa
airdata_all_wt_all <- subset(airdata_all_wt,type=="All")
wilcox.test(delta_pt~trt,data=airdata_all_wt_all)
airdata_all_wt_art <- subset(airdata_all_wt,type=="a")
wilcox.test(delta_pt~trt,data=airdata_all_wt_art)
airdata_all_wt_cap <- subset(airdata_all_wt,type=="c")
wilcox.test(delta_pt~trt,data=airdata_all_wt_cap)
airdata_all_wt_ven <- subset(airdata_all_wt,type=="v")
wilcox.test(delta_pt~trt,data=airdata_all_wt_ven)
label1.df <- data.frame(type = factor(c("All","v")),
                        delta_pt = c(1.5,2.1))
data1 <- data.frame(x = c(0.77,0.77,1.23,1.23), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data4 <- data.frame(x = c(3.77,3.77,4.23,4.23), y = c(2,2.1,2.1,2),
                    type=NA, delta_pt=NA)
aa<-aa + geom_path(data = data1, aes(x = x, y = y),size=1)+
  geom_path(data = data4, aes(x = x, y = y),size=1)+
  geom_text(data=label1.df,label="**",size=15)
aa

# 4 bars in each group
aa<-ggplot(data_summ, aes(x=type, y=delta_pt)) +
  geom_bar(aes(fill=gntrt),position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(fill=gntrt,ymin=delta_pt, ymax=delta_pt+ci),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Vessel Type") +
  ylab("??TTP (s)") + ylim(0,2.5)+
  scale_fill_manual(name="Genotype&Treatmet",values=c("greenyellow","green4","plum4","red4"),labels=c("nTg&nTx","Tg&nTx","nTg&Tx","Tg&Tx"))+
  scale_x_discrete("Vessel Type",limits=c("All","a","c","v"),
                   labels=c("All","Arterioles","Capillaries","Venules"))+
  theme(axis.text=element_text(size=20,colour="black"),
        text=element_text(size=20),
        legend.text=element_text(size=15))
label1.df <- data.frame(type = factor(c("All", "a")),
                        delta_pt = c(1.5,1.6))
label2.df <- data.frame(type = factor(c("v")),
                        delta_pt = c(2.1))
data1 <- data.frame(x = c(0.77,0.77,1.23,1.23), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(1.77,1.77,2.23,2.23), y = c(1.5, 1.6, 1.6, 1.5),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(3.77,3.77,4.23,4.23), y = c(2,2.1,2.1,2),
                    type=NA, delta_pt=NA)
aa<-aa + geom_path(data = data1, aes(x = x, y = y),size=1)+
  geom_path(data = data2, aes(x = x, y = y),size=1)+
  geom_path(data = data3, aes(x = x, y = y),size=1)+
  geom_text(data=label1.df,label="**",size=15)+
  geom_text(data=label2.df,label="*",size=15)
aa
dev.copy(jpeg,"??TTP_all.jpeg")
dev.off()
airdata_all_all_all <- subset(airdata_all_all,type=="All")
wilcox.test(delta_pt~gn,data=airdata_all_all_all)
airdata_all_all_art <- subset(airdata_all_all,type=="a")
wilcox.test(delta_pt~gn,data=airdata_all_all_art)
airdata_all_all_cap <- subset(airdata_all_all,type=="c")
wilcox.test(delta_pt~gn,data=airdata_all_all_cap)
airdata_all_all_ven <- subset(airdata_all_all,type=="v")
wilcox.test(delta_pt~gn,data=airdata_all_all_ven)





cc<-ggplot(data = airdata_ntx,aes(x=type,y=delta_pt))+geom_bar(stat = "identity")
cc+labs(y=expression(paste(Delta,"TTP (s)")))+
  labs(x="Vessel Type")+
  theme(axis.text=element_text(size=13),
        text=element_text(size=18),
        legend.text=element_text(size=15))+
  ggtitle(paste("All types of vessles"))+
  scale_x_discrete(labels=c("nTg","Tg"))
dev.copy(jpeg,"box_byGENO_ALL.jpeg")
dev.off()

data_ntx = subset(data_matched,trt=="ntx"&type=="a")
airdata_ntx <- subset(data_ntx,state=="air")
cdata_ntx <- subset(data_ntx,state=="co2")
delta_ntx_pt <- cdata_ntx$pt - airdata_ntx$pt
airdata_ntx$delta_pt <- delta_ntx_pt
airdata_ntx <- airdata_ntx[order(rev(airdata_ntx$gn),airdata_ntx$id,airdata_ntx$vessel_number),] #
airdata_ntx <- subset(airdata_ntx,id!=74)
airdata_ntx$gn <- factor(airdata_ntx$gn,c("wt","tg"))
cc<-ggplot(data = airdata_ntx,aes(x=gn,y=delta_pt,colour=gn))+geom_boxplot()
cc+labs(y=expression(paste(Delta,"TTP (s)")))+
  labs(x="Genotype")+
  theme(axis.text=element_text(size=13),
        text=element_text(size=18),
        legend.text=element_text(size=15))+
  ggtitle(paste("arterioles"))
dev.copy(jpeg,"box_byGENO_art.jpeg")
dev.off()

data_ntx = subset(data_matched,trt=="ntx"&type=="c")
airdata_ntx <- subset(data_ntx,state=="air")
cdata_ntx <- subset(data_ntx,state=="co2")
delta_ntx_pt <- cdata_ntx$pt - airdata_ntx$pt
airdata_ntx$delta_pt <- delta_ntx_pt
airdata_ntx <- airdata_ntx[order(rev(airdata_ntx$gn),airdata_ntx$id,airdata_ntx$vessel_number),] #
airdata_ntx <- subset(airdata_ntx,id!=74)
airdata_ntx$gn <- factor(airdata_ntx$gn,c("wt","tg"))
cc<-ggplot(data = airdata_ntx,aes(x=gn,y=delta_pt,colour=gn))+geom_boxplot()
cc+labs(y=expression(paste(Delta,"TTP (s)")))+
  labs(x="Genotype")+
  theme(axis.text=element_text(size=13),
        text=element_text(size=18),
        legend.text=element_text(size=15))+
  ggtitle(paste("capillaries"))
dev.copy(jpeg,"box_byGENO_cap.jpeg")
dev.off()

aa <-data.frame(airdata$pt,cdata$pt)
ggplot(aa,aes(x=airdata.pt,y=cdata.pt))+geom_point()+stat_smooth(method = 'lm',col="red")

## Linear mixed effects
# effects of genotype on resting TTP
data_ntx = subset(data_matched,trt=="ntx")
summary(lme(pt ~ state, random = list(~1|id), na.action = na.omit, data_ntx))
# INV. PROP.      Value Std.Error  DF   t-value p-value
# stateco2    -0.978116 0.1910925  68 -5.118547       0
data_ntx_flow = subset(data_matched,trt=="ntx" & slope!=0)
summary(lme(slope ~ gn, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,state=="air")))
# NO EFFECT       Value  Std.Error  DF  t-value p-value
# stateco2   0.05968779 0.05033579  13 1.185792  0.2569 !! NO genotypic Effect on slope during air !!
summary(lme(slope ~ gn, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,state=="co2")))
# NO EFFECT        Value  Std.Error  DF  t-value p-value
# stateco2    0.00943506 0.04854270  13 0.194366  0.8489 !! No genotypic Effect on slope during air !!
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,gn=="wt"&type=="a")))
# NO EFFECT         Value  Std.Error  DF   t-value p-value
# stateco2    -0.02679036 0.02964282  60 -0.903772  0.3697 !! No change on slope from air to co2 !!

summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,gn=="tg"&type=="a")))



summary(lme(AUC ~ gn, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,state=="air")))
# NO EFFECT       Value  Std.Error  DF  t-value p-value
# gnwt         271090.2   308196.7  13 0.879601   0.395 !! NO genotypic Effect on Area Under Curve during air !!
summary(lme(AUC ~ gn, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,state=="co2")))
# NO EFFECT        Value  Std.Error  DF   t-value p-value
# stateco2     -237604.2   180120.0  13 -1.319144  0.2099 !! No genotypic Effect on Area Under Curve during air !!
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,gn=="wt")))
# INV. PROP?      Value Std.Error  DF   t-value p-value
# stateco2      -382012  160718.7  60 -2.376898  0.0207
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,gn=="tg")))


# for RAT 150 only

summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==150)))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==150)))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==150)))
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==150 & type=="a")))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==150 & type=="a")))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==150 & type=="a")))
# for RAT 151 only
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==151)))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==151)))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==151)))
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==151 & type=="a")))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==151 & type=="a")))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==151 & type=="a")))
# for RAT 120 only
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==120)))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==120)))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==120)))
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==120 & type=="a")))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==120 & type=="a")))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==120 & type=="a")))
# for RAT 119 only
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==119)))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==119)))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==119)))
summary(lme(pt ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx,id==119 & type=="a")))
summary(lme(slope ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==119 & type=="a")))
summary(lme(AUC ~ state, random = list(~1|id, ~1|bolus_number), na.action = na.omit, subset(data_ntx_flow,id==119 & type=="a")))



### REGRESSIONS
# CI for Orthogonal Least Square fitting of a line: http://stackoverflow.com/questions/31057192/resampling-not-producing-expected-result-of-principal-component-analysis
## nTG nTX - matching air&co2
air_not_c <- airdata_wt_ntx[airdata_wt_ntx$type!="c",]
c_not_c <- cdata_wt_ntx[cdata_wt_ntx$type!="c",]
data_c_air <- subset(airdata_wt_ntx, type=="c")
data_c_co2 <- subset(cdata_wt_ntx,  type=="c")
keys <- c("vessel_number", "depth")
data_c_air_match <- merge(data_c_air, data_c_co2,by=keys)
ves <- data_c_air_match$vessel_number
data_c_air_mtchd <- subset(data_c_air, vessel_number %in% ves)
data_c_co2_mtchd <- subset(data_c_co2, vessel_number %in% ves)
data_c_mtchd <- rbind(data_c_air_mtchd,data_c_co2_mtchd)
data.frame(table(data_c_mtchd$state))   # number of air and co2 should be the same
airdata_wt_ntx_mtchd <- rbind(air_not_c,data_c_air_mtchd)
cdata_wt_ntx_mtchd <- rbind(c_not_c,data_c_co2_mtchd)
## TG nTX - matching air&co2
air_not_c <- airdata_tg_ntx[airdata_tg_ntx$type!="c",]
c_not_c <- cdata_tg_ntx[cdata_tg_ntx$type!="c",]
data_c_air <- subset(airdata_tg_ntx, type=="c")
data_c_co2 <- subset(cdata_tg_ntx,  type=="c")
keys <- c("vessel_number", "depth")
data_c_air_match <- merge(data_c_air, data_c_co2,by=keys)
ves <- data_c_air_match$vessel_number
data_c_air_mtchd <- subset(data_c_air, vessel_number %in% ves)
data_c_co2_mtchd <- subset(data_c_co2, vessel_number %in% ves)
data_c_mtchd <- rbind(data_c_air_mtchd,data_c_co2_mtchd)
data.frame(table(data_c_mtchd$state))   # number of air and co2 should be the same
airdata_tg_ntx_mtchd <- rbind(air_not_c,data_c_air_mtchd)
cdata_tg_ntx_mtchd <- rbind(c_not_c,data_c_co2_mtchd)

# Deming regression
f <- function(formula,data,SDR=1,...){
  M <- model.frame(formula, data)
  d <- Deming(x =M[,2],y =M[,1], sdr=SDR)[1:2]
  class(d) <- "Deming"
  d  
}
# an s3 method for predictdf (called within stat_smooth)
predictdf.Deming <- function(model, xseq, se, level) {
  pred <- model %*% t(cbind(1, xseq) )
  data.frame(x = xseq, y = c(pred))
}

# histogram showing nTg vs Tg to explain RFC
wt_ntx_all <-rbind(airdata_wt_ntx_mtchd,cdata_wt_ntx_mtchd)
aa = ggplot(wt_ntx_all,aes(pt,fill=state))+geom_histogram(alpha=0.5,aes(y=..density..),position='identity',
                                                          binwidth=0.5)
aa = ggplot(wt_ntx_all,aes(pt,fill=state))+geom_density(alpha=0.2)
aa
wt_ntx_art <-subset(wt_ntx_all,type=="a")
aa = ggplot(wt_ntx_art,aes(pt,fill=state))+geom_histogram(alpha=0.5,aes(y=..density..),position='identity',
                                                          binwidth=0.2)
aa

# all_regression wt ntx
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
colnames(wt_ntx_ac) <- c("air","co2")
wt_ntx_pca <- deming(wt_ntx_ac[,2]~wt_ntx_ac[,1])
wt_ntx_flow_change <- 1/wt_ntx_pca$coefficients[2]
wt_ntx_ci_int <- (wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1])/2
wt_ntx_flow_change_ci <- (wt_ntx_ci_int*wt_ntx_flow_change)^2
wt_ntx_slp_se <- sqrt(diag(wt_ntx_pca$variance))[2]
line_cu <- list(slope=wt_ntx_pca$ci[2,2],intercept=wt_ntx_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_pca$ci[2,1],intercept=wt_ntx_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_ac$cu <- wt_ntx_pca$ci[2,2]*wt_ntx_ac$air+wt_ntx_pca$ci[1,1]
wt_ntx_ac$cl <- wt_ntx_pca$ci[2,1]*wt_ntx_ac$air+wt_ntx_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_pca$ci[1,2]-wt_ntx_pca$ci[1,1])/(wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1]),
               y=(wt_ntx_pca$ci[2,1]*wt_ntx_pca$ci[1,1]-wt_ntx_pca$ci[2,2]*wt_ntx_pca$ci[1,2])/(wt_ntx_pca$ci[2,1]-wt_ntx_pca$ci[2,2])))
wt_ntx_plot <- ggplotRegressionPCA(wt_ntx_ac,5.5,12)
wt_ntx_plot
aa <- wt_ntx_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_nTG_nTX_ALL.jpeg")
dev.off()
wt_ntx_flow_change_tab = paste(round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_change_ci,digits=3))
# all_regression tg ntx
tg_ntx_ac <-data.frame(airdata_tg_ntx_mtchd$pt,cdata_tg_ntx_mtchd$pt)
colnames(tg_ntx_ac) <- c("air","co2")
tg_ntx_pca <- deming(tg_ntx_ac[,2]~tg_ntx_ac[,1])
tg_ntx_flow_change <- 1/tg_ntx_pca$coefficients[2]
tg_ntx_ci_int <- (tg_ntx_pca$ci[2,2]-tg_ntx_pca$ci[2,1])/2
tg_ntx_flow_change_ci <- (tg_ntx_ci_int*tg_ntx_flow_change)^2
tg_ntx_slp_se <- sqrt(diag(tg_ntx_pca$variance))[2]
line_cu <- list(slope=tg_ntx_pca$ci[2,2],intercept=tg_ntx_pca$ci[1,1])
line_cl <- list(slope=tg_ntx_pca$ci[2,1],intercept=tg_ntx_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,9,by=0.1))
colnames(lines) <- "x"
tg_ntx_ac$cu <- tg_ntx_pca$ci[2,2]*tg_ntx_ac$air+tg_ntx_pca$ci[1,1]
tg_ntx_ac$cl <- tg_ntx_pca$ci[2,1]*tg_ntx_ac$air+tg_ntx_pca$ci[1,2]
cu_coord <- tg_ntx_pca$ci[2,2]*x_coord+tg_ntx_pca$ci[1,1]
cl_coord <- tg_ntx_pca$ci[2,1]*x_coord+tg_ntx_pca$ci[1,2]
intersect <- as.data.frame(c(x=(tg_ntx_pca$ci[1,2]-tg_ntx_pca$ci[1,1])/(tg_ntx_pca$ci[2,2]-tg_ntx_pca$ci[2,1]),
                             y=(tg_ntx_pca$ci[2,1]*tg_ntx_pca$ci[1,1]-tg_ntx_pca$ci[2,2]*tg_ntx_pca$ci[1,2])/(tg_ntx_pca$ci[2,1]-tg_ntx_pca$ci[2,2])))
tg_ntx_plot <- ggplotRegressionPCA(tg_ntx_ac,5.5,12)
aa <- tg_ntx_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_TG_nTX_ALL.jpeg")
dev.off()
tg_ntx_flow_change_tab = paste(round(tg_ntx_flow_change,digits=3)," +/- ",round(tg_ntx_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_flow_change = (tg_ntx_pca$coefficient[2]-wt_ntx_pca$coefficient[2])/
  (sqrt(tg_ntx_slp_se^2+wt_ntx_slp_se^2))
p_ntx_flow_change=2*pnorm(-abs(z_ntx_flow_change))


# Arterioles_regression wt ntx
airdata_wt_ntx_art_mtchd<-subset(airdata_wt_ntx_mtchd,type=="a")
cdata_wt_ntx_art_mtchd<-subset(cdata_wt_ntx_mtchd,type=="a")
wt_ntx_art_ac <-data.frame(airdata_wt_ntx_art_mtchd$pt,cdata_wt_ntx_art_mtchd$pt)
colnames(wt_ntx_art_ac) <- c("air","co2")
wt_ntx_art_pca <- deming(wt_ntx_art_ac[,2]~wt_ntx_art_ac[,1])
wt_ntx_art_flow_change <- 1/wt_ntx_art_pca$coefficients[2]
wt_ntx_art_ci_int <- (wt_ntx_art_pca$ci[2,2]-wt_ntx_art_pca$ci[2,1])/2
wt_ntx_art_flow_change_ci <- (wt_ntx_art_ci_int*wt_ntx_art_flow_change)^2
wt_ntx_art_slp_se <- sqrt(diag(wt_ntx_art_pca$variance))[2]
line_cu <- list(slope=wt_ntx_art_pca$ci[2,2],intercept=wt_ntx_art_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_art_pca$ci[2,1],intercept=wt_ntx_art_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_art_ac$cu <- wt_ntx_art_pca$ci[2,2]*wt_ntx_art_ac$air+wt_ntx_art_pca$ci[1,1]
wt_ntx_art_ac$cl <- wt_ntx_art_pca$ci[2,1]*wt_ntx_art_ac$air+wt_ntx_art_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_art_pca$ci[1,2]-wt_ntx_art_pca$ci[1,1])/(wt_ntx_art_pca$ci[2,2]-wt_ntx_art_pca$ci[2,1]),
                             y=(wt_ntx_art_pca$ci[2,1]*wt_ntx_art_pca$ci[1,1]-wt_ntx_art_pca$ci[2,2]*wt_ntx_art_pca$ci[1,2])/(wt_ntx_art_pca$ci[2,1]-wt_ntx_art_pca$ci[2,2])))
wt_ntx_art_plot <- ggplotRegressionPCA(wt_ntx_art_ac,5.5,12)
aa <- wt_ntx_art_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Arterioles")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_nTG_nTX_Arterioles.jpeg")
dev.off()
wt_ntx_art_flow_change_tab = paste(round(wt_ntx_art_flow_change,digits=3)," +/- ",round(wt_ntx_art_flow_change_ci,digits=3))
# Arterioles_regression tg ntx
airdata_tg_ntx_art_mtchd<-subset(airdata_tg_ntx_mtchd,type=="a")
cdata_tg_ntx_art_mtchd<-subset(cdata_tg_ntx_mtchd,type=="a")
tg_ntx_art_ac <-data.frame(airdata_tg_ntx_art_mtchd$pt,cdata_tg_ntx_art_mtchd$pt)
colnames(tg_ntx_art_ac) <- c("air","co2")
tg_ntx_art_pca <- deming(tg_ntx_art_ac[,2]~tg_ntx_art_ac[,1])
tg_ntx_art_flow_change <- 1/tg_ntx_art_pca$coefficients[2]
tg_ntx_art_ci_int <- (tg_ntx_art_pca$ci[2,2]-tg_ntx_art_pca$ci[2,1])/2
tg_ntx_art_flow_change_ci <- (tg_ntx_art_ci_int*tg_ntx_art_flow_change)^2
tg_ntx_art_slp_se <- sqrt(diag(tg_ntx_art_pca$variance))[2]
line_cu <- list(slope=tg_ntx_art_pca$ci[2,2],intercept=tg_ntx_art_pca$ci[1,1])
line_cl <- list(slope=tg_ntx_art_pca$ci[2,1],intercept=tg_ntx_art_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
tg_ntx_art_ac$cu <- tg_ntx_art_pca$ci[2,2]*tg_ntx_art_ac$air+tg_ntx_art_pca$ci[1,1]
tg_ntx_art_ac$cl <- tg_ntx_art_pca$ci[2,1]*tg_ntx_art_ac$air+tg_ntx_art_pca$ci[1,2]
intersect <- as.data.frame(c(x=(tg_ntx_art_pca$ci[1,2]-tg_ntx_art_pca$ci[1,1])/(tg_ntx_art_pca$ci[2,2]-tg_ntx_art_pca$ci[2,1]),
                             y=(tg_ntx_art_pca$ci[2,1]*tg_ntx_art_pca$ci[1,1]-tg_ntx_art_pca$ci[2,2]*tg_ntx_art_pca$ci[1,2])/(tg_ntx_art_pca$ci[2,1]-tg_ntx_art_pca$ci[2,2])))
tg_ntx_art_plot <- ggplotRegressionPCA(tg_ntx_art_ac,5.5,12)
aa <- tg_ntx_art_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Arterioles")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_TG_nTX_Arterioles.jpeg")
dev.off()
tg_ntx_art_flow_change_tab = paste(round(tg_ntx_art_flow_change,digits=3)," +/- ",round(tg_ntx_art_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_art_flow_change = (tg_ntx_art_pca$coefficient[2]-wt_ntx_art_pca$coefficient[2])/
  (sqrt(tg_ntx_art_slp_se^2+wt_ntx_art_slp_se^2))
p_ntx_art_flow_change=2*pnorm(-abs(z_ntx_art_flow_change))


# Capillaries_regression wt ntx
airdata_wt_ntx_cap_mtchd<-subset(airdata_wt_ntx_mtchd,type=="c")
cdata_wt_ntx_cap_mtchd<-subset(cdata_wt_ntx_mtchd,type=="c")
wt_ntx_cap_ac <-data.frame(airdata_wt_ntx_cap_mtchd$pt,cdata_wt_ntx_cap_mtchd$pt)
colnames(wt_ntx_cap_ac) <- c("air","co2")
wt_ntx_cap_pca <- deming(wt_ntx_cap_ac[,2]~wt_ntx_cap_ac[,1])
wt_ntx_cap_flow_change <- 1/wt_ntx_cap_pca$coefficients[2]
wt_ntx_cap_ci_int <- (wt_ntx_cap_pca$ci[2,2]-wt_ntx_cap_pca$ci[2,1])/2
wt_ntx_cap_flow_change_ci <- (wt_ntx_cap_ci_int*wt_ntx_cap_flow_change)^2
wt_ntx_cap_slp_se <- sqrt(diag(wt_ntx_cap_pca$variance))[2]
line_cu <- list(slope=wt_ntx_cap_pca$ci[2,2],intercept=wt_ntx_cap_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_cap_pca$ci[2,1],intercept=wt_ntx_cap_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_cap_ac$cu <- wt_ntx_cap_pca$ci[2,2]*wt_ntx_cap_ac$air+wt_ntx_cap_pca$ci[1,1]
wt_ntx_cap_ac$cl <- wt_ntx_cap_pca$ci[2,1]*wt_ntx_cap_ac$air+wt_ntx_cap_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_cap_pca$ci[1,2]-wt_ntx_cap_pca$ci[1,1])/(wt_ntx_cap_pca$ci[2,2]-wt_ntx_cap_pca$ci[2,1]),
                             y=(wt_ntx_cap_pca$ci[2,1]*wt_ntx_cap_pca$ci[1,1]-wt_ntx_cap_pca$ci[2,2]*wt_ntx_cap_pca$ci[1,2])/(wt_ntx_cap_pca$ci[2,1]-wt_ntx_cap_pca$ci[2,2])))
wt_ntx_cap_plot <- ggplotRegressionPCA(wt_ntx_cap_ac,5.5,12)
aa <- wt_ntx_cap_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Capillaries")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_nTG_nTX_Capillaries.jpeg")
dev.off()
wt_ntx_cap_flow_change_tab = paste(round(wt_ntx_cap_flow_change,digits=3)," +/- ",round(wt_ntx_cap_flow_change_ci,digits=3))
# Capillaries_regression tg ntx
airdata_tg_ntx_cap_mtchd<-subset(airdata_tg_ntx_mtchd,type=="c")
cdata_tg_ntx_cap_mtchd<-subset(cdata_tg_ntx_mtchd,type=="c")
tg_ntx_cap_ac <-data.frame(airdata_tg_ntx_cap_mtchd$pt,cdata_tg_ntx_cap_mtchd$pt)
colnames(tg_ntx_cap_ac) <- c("air","co2")
tg_ntx_cap_pca <- deming(tg_ntx_cap_ac[,2]~tg_ntx_cap_ac[,1])
tg_ntx_cap_flow_change <- 1/tg_ntx_cap_pca$coefficients[2]
tg_ntx_cap_ci_int <- (tg_ntx_cap_pca$ci[2,2]-tg_ntx_cap_pca$ci[2,1])/2
tg_ntx_cap_flow_change_ci <- (tg_ntx_cap_ci_int*tg_ntx_cap_flow_change)^2
tg_ntx_cap_slp_se <- sqrt(diag(tg_ntx_cap_pca$variance))[2]
line_cu <- list(slope=tg_ntx_cap_pca$ci[2,2],intercept=tg_ntx_cap_pca$ci[1,1])
line_cl <- list(slope=tg_ntx_cap_pca$ci[2,1],intercept=tg_ntx_cap_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
tg_ntx_cap_ac$cu <- tg_ntx_cap_pca$ci[2,2]*tg_ntx_cap_ac$air+tg_ntx_cap_pca$ci[1,1]
tg_ntx_cap_ac$cl <- tg_ntx_cap_pca$ci[2,1]*tg_ntx_cap_ac$air+tg_ntx_cap_pca$ci[1,2]
intersect <- as.data.frame(c(x=(tg_ntx_cap_pca$ci[1,2]-tg_ntx_cap_pca$ci[1,1])/(tg_ntx_cap_pca$ci[2,2]-tg_ntx_cap_pca$ci[2,1]),
                             y=(tg_ntx_cap_pca$ci[2,1]*tg_ntx_cap_pca$ci[1,1]-tg_ntx_cap_pca$ci[2,2]*tg_ntx_cap_pca$ci[1,2])/(tg_ntx_cap_pca$ci[2,1]-tg_ntx_cap_pca$ci[2,2])))
tg_ntx_cap_plot <- ggplotRegressionPCA(tg_ntx_cap_ac,5.5,12)
aa <- tg_ntx_cap_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Capillaries")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_TG_nTX_Capillaries.jpeg")
dev.off()
tg_ntx_cap_flow_change_tab = paste(round(tg_ntx_cap_flow_change,digits=3)," +/- ",round(tg_ntx_cap_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_cap_flow_change = (tg_ntx_cap_pca$coefficient[2]-wt_ntx_cap_pca$coefficient[2])/
  (sqrt(tg_ntx_cap_slp_se^2+wt_ntx_cap_slp_se^2))
p_ntx_cap_flow_change=2*pnorm(-abs(z_ntx_cap_flow_change))


# Venules_regression wt ntx
airdata_wt_ntx_ven_mtchd<-subset(airdata_wt_ntx_mtchd,type=="v")
cdata_wt_ntx_ven_mtchd<-subset(cdata_wt_ntx_mtchd,type=="v")
wt_ntx_ven_ac <-data.frame(airdata_wt_ntx_ven_mtchd$pt,cdata_wt_ntx_ven_mtchd$pt)
colnames(wt_ntx_ven_ac) <- c("air","co2")
wt_ntx_ven_pca <- deming(wt_ntx_ven_ac[,2]~wt_ntx_ven_ac[,1])
wt_ntx_ven_flow_change <- 1/wt_ntx_ven_pca$coefficients[2]
wt_ntx_ven_ci_int <- (wt_ntx_ven_pca$ci[2,2]-wt_ntx_ven_pca$ci[2,1])/2
wt_ntx_ven_flow_change_ci <- (wt_ntx_ven_ci_int*wt_ntx_ven_flow_change)^2
wt_ntx_ven_slp_se <- sqrt(diag(wt_ntx_ven_pca$variance))[2]
line_cu <- list(slope=wt_ntx_ven_pca$ci[2,2],intercept=wt_ntx_ven_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_ven_pca$ci[2,1],intercept=wt_ntx_ven_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_ven_ac$cu <- wt_ntx_ven_pca$ci[2,2]*wt_ntx_ven_ac$air+wt_ntx_ven_pca$ci[1,1]
wt_ntx_ven_ac$cl <- wt_ntx_ven_pca$ci[2,1]*wt_ntx_ven_ac$air+wt_ntx_ven_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_ven_pca$ci[1,2]-wt_ntx_ven_pca$ci[1,1])/(wt_ntx_ven_pca$ci[2,2]-wt_ntx_ven_pca$ci[2,1]),
                             y=(wt_ntx_ven_pca$ci[2,1]*wt_ntx_ven_pca$ci[1,1]-wt_ntx_ven_pca$ci[2,2]*wt_ntx_ven_pca$ci[1,2])/(wt_ntx_ven_pca$ci[2,1]-wt_ntx_ven_pca$ci[2,2])))
wt_ntx_ven_plot <- ggplotRegressionPCA(wt_ntx_ven_ac,5.5,12)
aa <- wt_ntx_ven_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Venules")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_nTG_nTX_Venules.jpeg")
dev.off()
wt_ntx_ven_flow_change_tab = paste(round(wt_ntx_ven_flow_change,digits=3)," +/- ",round(wt_ntx_ven_flow_change_ci,digits=3))
# Venules_regression tg ntx
airdata_tg_ntx_ven_mtchd<-subset(airdata_tg_ntx_mtchd,type=="v")
cdata_tg_ntx_ven_mtchd<-subset(cdata_tg_ntx_mtchd,type=="v")
tg_ntx_ven_ac <-data.frame(airdata_tg_ntx_ven_mtchd$pt,cdata_tg_ntx_ven_mtchd$pt)
colnames(tg_ntx_ven_ac) <- c("air","co2")
tg_ntx_ven_pca <- deming(tg_ntx_ven_ac[,2]~tg_ntx_ven_ac[,1])
tg_ntx_ven_flow_change <- 1/tg_ntx_ven_pca$coefficients[2]
tg_ntx_ven_ci_int <- (tg_ntx_ven_pca$ci[2,2]-tg_ntx_ven_pca$ci[2,1])/2
tg_ntx_ven_flow_change_ci <- (tg_ntx_ven_ci_int*tg_ntx_ven_flow_change)^2
tg_ntx_ven_slp_se <- sqrt(diag(tg_ntx_ven_pca$variance))[2]
line_cu <- list(slope=tg_ntx_ven_pca$ci[2,2],intercept=tg_ntx_ven_pca$ci[1,1])
line_cl <- list(slope=tg_ntx_ven_pca$ci[2,1],intercept=tg_ntx_ven_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
tg_ntx_ven_ac$cu <- tg_ntx_ven_pca$ci[2,2]*tg_ntx_ven_ac$air+tg_ntx_ven_pca$ci[1,1]
tg_ntx_ven_ac$cl <- tg_ntx_ven_pca$ci[2,1]*tg_ntx_ven_ac$air+tg_ntx_ven_pca$ci[1,2]
intersect <- as.data.frame(c(x=(tg_ntx_ven_pca$ci[1,2]-tg_ntx_ven_pca$ci[1,1])/(tg_ntx_ven_pca$ci[2,2]-tg_ntx_ven_pca$ci[2,1]),
                             y=(tg_ntx_ven_pca$ci[2,1]*tg_ntx_ven_pca$ci[1,1]-tg_ntx_ven_pca$ci[2,2]*tg_ntx_ven_pca$ci[1,2])/(tg_ntx_ven_pca$ci[2,1]-tg_ntx_ven_pca$ci[2,2])))
tg_ntx_ven_plot <- ggplotRegressionPCA(tg_ntx_ven_ac,5.5,12)
aa <- tg_ntx_ven_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Venules")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_AUC_TG_nTX_Venules.jpeg")
dev.off()
tg_ntx_ven_flow_change_tab = paste(round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_ven_flow_change = (tg_ntx_ven_pca$coefficient[2]-wt_ntx_ven_pca$coefficient[2])/
  (sqrt(tg_ntx_ven_slp_se^2+wt_ntx_ven_slp_se^2))
p_ntx_ven_flow_change=2*pnorm(-abs(z_ntx_ven_flow_change))

wt_ntx_flow_tabs <- data.frame(c(wt_ntx_flow_change_tab,wt_ntx_art_flow_change_tab,wt_ntx_cap_flow_change_tab,wt_ntx_ven_flow_change_tab))
tg_ntx_flow_tabs <- data.frame(c(tg_ntx_flow_change_tab,tg_ntx_art_flow_change_tab,tg_ntx_cap_flow_change_tab,tg_ntx_ven_flow_change_tab))
p_values <- data.frame(c(p_ntx_flow_change,p_ntx_art_flow_change,p_ntx_cap_flow_change,p_ntx_ven_flow_change))
flow_change_tabs <- cbind(wt_ntx_flow_tabs,tg_ntx_flow_tabs,p_values)
colnames(flow_change_tabs) <- c("nTg Flow Change","Tg Flow Change","p-value")
rownames(flow_change_tabs) <- c("ALL","Arterioles","Capillaries","Venules")
flow_change_tabs

# /Ordinary regression (error in y only) -----------------------------------


slope=wt_ntx_pca$ci[2,2],intercept=wt_ntx_pca$ci[1,1]
geom_segment(x=min(x_coord),y=min(cu_coord),xend=max(x_coord),yend=max(cu_coord),
             colour='darkgrey')+
  geom_segment(x=min(x_coord),y=min(cl_coord),xend=max(x_coord),yend=max(cl_coord),
               colour='darkgrey')+

wt_ntx_pca <- prcomp(wt_ntx_ac)$rotation
wt_ntx_pca_slp <- wt_ntx_pca[2,1]/wt_ntx_pca[1,1]
set.seed(1)
mydf <- wt_ntx_ac[1:2]
times <- 1000
ll <- vector(mode = "list", length = times)
for (i in seq_len(times)) {
  tempdf  <- mydf[sample(nrow(mydf), replace = TRUE), ]
  ll[[i]] <- abs(prcomp(tempdf)$rotation) ## NOTE: abs(...)
}
wt_ntx_pca_cl <- data.frame(apply(simplify2array(ll), 1:2, quantile, probs = 0.025))
wt_ntx_pca_slp_cl <- wt_ntx_pca_cl[2,1]/wt_ntx_pca_cl[1,1]
wt_ntx_pca_cu <- data.frame(apply(simplify2array(ll), 1:2, quantile, probs = 0.975))
wt_ntx_pca_slp_cu <- wt_ntx_pca_cu[2,1]/wt_ntx_pca_cu[1,1]

wt_ntx_linfit <- lm(cdata_wt_ntx_mtchd.pt~airdata_wt_ntx_mtchd.pt,data=wt_ntx_ac)
wt_ntx_ci <- confint(wt_ntx_linfit, 'airdata_wt_ntx_mtchd.pt', level=0.95)
wt_ntx_ci_int <- (wt_ntx_ci[2]-wt_ntx_ci[1])/2
wt_ntx_gg <- ggplotRegression(wt_ntx_linfit)
wt_ntx_linfit_sm <-summary(wt_ntx_linfit)
wt_ntx_linfit_std <- wt_ntx_linfit_sm$coefficients[[4]]
wt_ntx_flow_change <- 1/wt_ntx_linfit$coefficients[2]
wt_ntx_flow_change_ci <- (wt_ntx_ci_int*wt_ntx_flow_change)^2
# wt_ntx_flow_change_std <- sqrt((wt_ntx_linfit_std*wt_ntx_flow_change/wt_ntx_linfit$coefficients[2])^2)
wt_ntx_gg_lab <- wt_ntx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg")
wt_ntx_gg_lab <- wt_ntx_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
wt_ntx_gg_lab <- wt_ntx_gg_lab+theme(axis.text=element_text(size=20),
                                          text=element_text(size=20))
wt_ntx_gg_lab
dev.copy(jpeg,"FC_nTG_nTX_rats.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX rats RFC = ",round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_change_ci,digits=3)," R2 = ",signif(summary(wt_ntx_linfit)$adj.r.squared, 3)))
wt_ntx_flow_change_tab = paste(round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_change_ci,digits=3))
# all_regression wt ntx
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
wt_ntx_linfit <- lm(cdata_wt_ntx_mtchd.pt~airdata_wt_ntx_mtchd.pt,data=wt_ntx_ac)
wt_ntx_ci <- confint(wt_ntx_linfit, 'airdata_wt_ntx_mtchd.pt', level=0.95)
wt_ntx_ci_int <- (wt_ntx_ci[2]-wt_ntx_ci[1])/2
wt_ntx_gg <- ggplotRegression(wt_ntx_linfit)
wt_ntx_linfit_sm <-summary(wt_ntx_linfit)
wt_ntx_linfit_std <- wt_ntx_linfit_sm$coefficients[[4]]
wt_ntx_flow_change <- 1/wt_ntx_linfit$coefficients[2]
wt_ntx_flow_change_ci <- (wt_ntx_ci_int*wt_ntx_flow_change)^2
# wt_ntx_flow_change_std <- sqrt((wt_ntx_linfit_std*wt_ntx_flow_change/wt_ntx_linfit$coefficients[2])^2)
wt_ntx_gg_lab <- wt_ntx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("wt")
wt_ntx_gg_lab <- wt_ntx_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
wt_ntx_gg_lab <- wt_ntx_gg_lab+theme(axis.text=element_text(size=20),
                                     text=element_text(size=20))
wt_ntx_gg_lab
dev.copy(jpeg,"FC_WT_nTX_rats.jpeg")
dev.off()
#+ggtitle(paste("TG nTX rats RFC = ",round(tg_ntx_flow_change,digits=2)," +/- ",round(tg_ntx_flow_change_ci,digits=3)," R2 = ",signif(summary(tg_ntx_linfit)$adj.r.squared, 3)))
wt_ntx_flow_change_tab = paste(round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_change_ci,digits=3))
# all_regression wt ntx
tg_ntx_ac <-data.frame(airdata_tg_ntx_mtchd$pt,cdata_tg_ntx_mtchd$pt)
tg_ntx_linfit <- lm(cdata_tg_ntx_mtchd.pt~airdata_tg_ntx_mtchd.pt,data=tg_ntx_ac)
tg_ntx_ci <- confint(tg_ntx_linfit, 'airdata_tg_ntx_mtchd.pt', level=0.95)
tg_ntx_ci_int <- (tg_ntx_ci[2]-tg_ntx_ci[1])/2
tg_ntx_gg <- ggplotRegression(tg_ntx_linfit)
tg_ntx_linfit_sm <-summary(tg_ntx_linfit)
tg_ntx_linfit_std <- tg_ntx_linfit_sm$coefficients[[4]]
tg_ntx_flow_change <- 1/tg_ntx_linfit$coefficients[2]
tg_ntx_flow_change_ci <- (tg_ntx_ci_int*tg_ntx_flow_change)^2
# tg_ntx_flow_change_std <- sqrt((tg_ntx_linfit_std*tg_ntx_flow_change/tg_ntx_linfit$coefficients[2])^2)
tg_ntx_gg_lab <- tg_ntx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg")
tg_ntx_gg_lab <- tg_ntx_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
tg_ntx_gg_lab <- tg_ntx_gg_lab+theme(axis.text=element_text(size=20),
                                     text=element_text(size=20))
tg_ntx_gg_lab
dev.copy(jpeg,"FC_TG_nTX_rats.jpeg")
dev.off()
#+ggtitle(paste("TG nTX rats RFC = ",round(tg_ntx_flow_change,digits=2)," +/- ",round(tg_ntx_flow_change_ci,digits=3)," R2 = ",signif(summary(tg_ntx_linfit)$adj.r.squared, 3)))
tg_ntx_flow_change_tab = paste(round(tg_ntx_flow_change,digits=3)," +/- ",round(tg_ntx_flow_change_ci,digits=3))
# z for hypothesis testing of two linear regression
# Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_flow_change = (tg_ntx_linfit$coefficients[2]-wt_ntx_linfit$coefficients[2])/
  (sqrt(tg_ntx_linfit_std^2+wt_ntx_linfit_std^2))
p_ntx_flow_change=2*pnorm(-abs(z_ntx_flow_change))

# arteriole_regression wt ntx
airdata_wt_ntx_art_mtchd<-subset(airdata_wt_ntx_mtchd,type=="a")
cdata_wt_ntx_art_mtchd<-subset(cdata_wt_ntx_mtchd,type=="a")
wt_ntx_art_ac <-data.frame(airdata_wt_ntx_art_mtchd$pt,cdata_wt_ntx_art_mtchd$pt)
wt_ntx_art_linfit <- lm(cdata_wt_ntx_art_mtchd.pt~airdata_wt_ntx_art_mtchd.pt,data=wt_ntx_art_ac)
wt_ntx_art_ci <- confint(wt_ntx_art_linfit, 'airdata_wt_ntx_art_mtchd.pt', level=0.95)
wt_ntx_art_ci_int <- (wt_ntx_art_ci[2]-wt_ntx_art_ci[1])/2
wt_ntx_art_gg <- ggplotRegression(wt_ntx_art_linfit)
wt_ntx_art_linfit_sm <-summary(wt_ntx_art_linfit)
wt_ntx_art_linfit_std <- wt_ntx_art_linfit_sm$coefficients[[4]]
wt_ntx_art_flow_change <- 1/wt_ntx_art_linfit$coefficients[2]
wt_ntx_art_flow_change_ci <- (wt_ntx_art_ci_int*wt_ntx_art_flow_change)^2
# wt_ntx_art_flow_change_std <- sqrt((wt_ntx_art_linfit_std*wt_ntx_art_flow_change/wt_ntx_art_linfit$coefficients[2])^2)
wt_ntx_art_gg_lab <- wt_ntx_art_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Arterioles")
wt_ntx_art_gg_lab <- wt_ntx_art_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
wt_ntx_art_gg_lab <- wt_ntx_art_gg_lab+theme(axis.text=element_text(size=20),
                                     text=element_text(size=20))
wt_ntx_art_gg_lab
dev.copy(jpeg,"FC_nTG_nTX_art.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX arterioles RFC = ",round(wt_ntx_art_flow_change,digits=3)," +/- ",round(wt_ntx_art_flow_change_ci,digits=3)," R2 = ",signif(summary(wt_ntx_art_linfit)$adj.r.squared, 3)))
wt_ntx_art_flow_change_tab = paste(round(wt_ntx_art_flow_change,digits=3)," +/- ",round(wt_ntx_art_flow_change_ci,digits=3))
# arterioles_regression tg ntx
airdata_tg_ntx_art_mtchd<-subset(airdata_tg_ntx_mtchd,type=="a")
cdata_tg_ntx_art_mtchd<-subset(cdata_tg_ntx_mtchd,type=="a")
tg_ntx_art_ac <-data.frame(airdata_tg_ntx_art_mtchd$pt,cdata_tg_ntx_art_mtchd$pt)
tg_ntx_art_linfit <- lm(cdata_tg_ntx_art_mtchd.pt~airdata_tg_ntx_art_mtchd.pt,data=tg_ntx_art_ac)
tg_ntx_art_ci <- confint(tg_ntx_art_linfit, 'airdata_tg_ntx_art_mtchd.pt', level=0.95)
tg_ntx_art_ci_int <- (tg_ntx_art_ci[2]-tg_ntx_art_ci[1])/2
tg_ntx_art_gg <- ggplotRegression(tg_ntx_art_linfit)
tg_ntx_art_linfit_sm <-summary(tg_ntx_art_linfit)
tg_ntx_art_linfit_std <- tg_ntx_art_linfit_sm$coefficients[[4]]
tg_ntx_art_flow_change <- 1/tg_ntx_art_linfit$coefficients[2]
tg_ntx_art_flow_change_ci <- (tg_ntx_art_ci_int*tg_ntx_art_flow_change)^2
# tg_ntx_art_flow_change_std <- sqrt((tg_ntx_art_linfit_std*tg_ntx_art_flow_change/tg_ntx_art_linfit$coefficients[2])^2)
tg_ntx_art_gg_lab <- tg_ntx_art_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Arterioles")
tg_ntx_art_gg_lab <- tg_ntx_art_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
tg_ntx_art_gg_lab <- tg_ntx_art_gg_lab+theme(axis.text=element_text(size=20),
                                             text=element_text(size=20))
tg_ntx_art_gg_lab
dev.copy(jpeg,"FC_TG_nTX_art.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX arterioles RFC = ",round(tg_ntx_art_flow_change,digits=3)," +/- ",round(tg_ntx_art_flow_change_ci,digits=3)," R2 = ",signif(summary(tg_ntx_art_linfit)$adj.r.squared, 3)))
tg_ntx_art_flow_change_tab = paste(round(tg_ntx_art_flow_change,digits=3)," +/- ",round(tg_ntx_art_flow_change_ci,digits=3))
# z for hypothesis testing of two linear regression
# Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_art_flow_change = (tg_ntx_art_linfit$coefficients[2]-wt_ntx_art_linfit$coefficients[2])/
  (sqrt(tg_ntx_art_linfit_std^2+wt_ntx_art_linfit_std^2))
p_ntx_art_flow_change=2*pnorm(-abs(z_ntx_art_flow_change))


# capillary_regression wt ntx
airdata_wt_ntx_cap_mtchd<-subset(airdata_wt_ntx_mtchd,type=="c")
cdata_wt_ntx_cap_mtchd<-subset(cdata_wt_ntx_mtchd,type=="c")
wt_ntx_cap_ac <-data.frame(airdata_wt_ntx_cap_mtchd$pt,cdata_wt_ntx_cap_mtchd$pt)
wt_ntx_cap_linfit <- lm(cdata_wt_ntx_cap_mtchd.pt~airdata_wt_ntx_cap_mtchd.pt,data=wt_ntx_cap_ac)
wt_ntx_cap_ci <- confint(wt_ntx_cap_linfit, 'airdata_wt_ntx_cap_mtchd.pt', level=0.95)
wt_ntx_cap_ci_int <- (wt_ntx_cap_ci[2]-wt_ntx_cap_ci[1])/2
wt_ntx_cap_gg <- ggplotRegression(wt_ntx_cap_linfit)
wt_ntx_cap_linfit_sm <-summary(wt_ntx_cap_linfit)
wt_ntx_cap_linfit_std <- wt_ntx_cap_linfit_sm$coefficients[[4]]
wt_ntx_cap_flow_change <- 1/wt_ntx_cap_linfit$coefficients[2]
wt_ntx_cap_flow_change_ci <- (wt_ntx_cap_ci_int*wt_ntx_cap_flow_change)^2
# wt_ntx_cap_flow_change_std <- sqrt((wt_ntx_cap_linfit_std*wt_ntx_cap_flow_change/wt_ntx_cap_linfit$coefficients[2])^2)
wt_ntx_cap_gg_lab <- wt_ntx_cap_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Capillaries")
wt_ntx_cap_gg_lab <- wt_ntx_cap_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
wt_ntx_cap_gg_lab <- wt_ntx_cap_gg_lab+theme(axis.text=element_text(size=20),
                                             text=element_text(size=20))
wt_ntx_cap_gg_lab
dev.copy(jpeg,"FC_nTG_nTX_cap.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX capillarys RFC = ",round(wt_ntx_cap_flow_change,digits=3)," +/- ",round(wt_ntx_cap_flow_change_ci,digits=3)," R2 = ",signif(summary(wt_ntx_cap_linfit)$adj.r.squared, 3)))
wt_ntx_cap_flow_change_tab = paste(round(wt_ntx_cap_flow_change,digits=3)," +/- ",round(wt_ntx_cap_flow_change_ci,digits=3))
# capillary_regression tg ntx
airdata_tg_ntx_cap_mtchd<-subset(airdata_tg_ntx_mtchd,type=="c")
cdata_tg_ntx_cap_mtchd<-subset(cdata_tg_ntx_mtchd,type=="c")
tg_ntx_cap_ac <-data.frame(airdata_tg_ntx_cap_mtchd$pt,cdata_tg_ntx_cap_mtchd$pt)
tg_ntx_cap_linfit <- lm(cdata_tg_ntx_cap_mtchd.pt~airdata_tg_ntx_cap_mtchd.pt,data=tg_ntx_cap_ac)
tg_ntx_cap_ci <- confint(tg_ntx_cap_linfit, 'airdata_tg_ntx_cap_mtchd.pt', level=0.95)
tg_ntx_cap_ci_int <- (tg_ntx_cap_ci[2]-tg_ntx_cap_ci[1])/2
tg_ntx_cap_gg <- ggplotRegression(tg_ntx_cap_linfit)
tg_ntx_cap_linfit_sm <-summary(tg_ntx_cap_linfit)
tg_ntx_cap_linfit_std <- tg_ntx_cap_linfit_sm$coefficients[[4]]
tg_ntx_cap_flow_change <- 1/tg_ntx_cap_linfit$coefficients[2]
tg_ntx_cap_flow_change_ci <- (tg_ntx_cap_ci_int*tg_ntx_cap_flow_change)^2
# tg_ntx_cap_flow_change_std <- sqrt((tg_ntx_cap_linfit_std*tg_ntx_cap_flow_change/tg_ntx_cap_linfit$coefficients[2])^2)
tg_ntx_cap_gg_lab <- tg_ntx_cap_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Capillaries")
tg_ntx_cap_gg_lab <- tg_ntx_cap_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
tg_ntx_cap_gg_lab <- tg_ntx_cap_gg_lab+theme(axis.text=element_text(size=20),
                                             text=element_text(size=20))
tg_ntx_cap_gg_lab
dev.copy(jpeg,"FC_TG_nTX_cap.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX arterioles RFC = ",round(tg_ntx_cap_flow_change,digits=3)," +/- ",round(tg_ntx_cap_flow_change_ci,digits=3)," R2 = ",signif(summary(tg_ntx_cap_linfit)$adj.r.squared, 3)))
tg_ntx_cap_flow_change_tab = paste(round(tg_ntx_cap_flow_change,digits=3)," +/- ",round(tg_ntx_cap_flow_change_ci,digits=3))
# z for hypothesis testing of two linear regression
# Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_cap_flow_change = (tg_ntx_cap_linfit$coefficients[2]-wt_ntx_cap_linfit$coefficients[2])/
  (sqrt(tg_ntx_cap_linfit_std^2+wt_ntx_cap_linfit_std^2))
p_ntx_cap_flow_change=2*pnorm(-abs(z_ntx_cap_flow_change))


# venule_regression wt ntx
airdata_wt_ntx_ven_mtchd<-subset(airdata_wt_ntx_mtchd,type=="v")
cdata_wt_ntx_ven_mtchd<-subset(cdata_wt_ntx_mtchd,type=="v")
wt_ntx_ven_ac <-data.frame(airdata_wt_ntx_ven_mtchd$pt,cdata_wt_ntx_ven_mtchd$pt)
wt_ntx_ven_linfit <- lm(cdata_wt_ntx_ven_mtchd.pt~airdata_wt_ntx_ven_mtchd.pt,data=wt_ntx_ven_ac)
wt_ntx_ven_ci <- confint(wt_ntx_ven_linfit, 'airdata_wt_ntx_ven_mtchd.pt', level=0.95)
wt_ntx_ven_ci_int <- (wt_ntx_ven_ci[2]-wt_ntx_ven_ci[1])/2
wt_ntx_ven_gg <- ggplotRegression(wt_ntx_ven_linfit)
wt_ntx_ven_linfit_sm <-summary(wt_ntx_ven_linfit)
wt_ntx_ven_linfit_std <- wt_ntx_ven_linfit_sm$coefficients[[4]]
wt_ntx_ven_flow_change <- 1/wt_ntx_ven_linfit$coefficients[2]
wt_ntx_ven_flow_change_ci <- (wt_ntx_ven_ci_int*wt_ntx_ven_flow_change)^2
# wt_ntx_ven_flow_change_std <- sqrt((wt_ntx_ven_linfit_std*wt_ntx_ven_flow_change/wt_ntx_ven_linfit$coefficients[2])^2)
wt_ntx_ven_gg_lab <- wt_ntx_ven_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - Venules")
wt_ntx_ven_gg_lab <- wt_ntx_ven_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
wt_ntx_ven_gg_lab <- wt_ntx_ven_gg_lab+theme(axis.text=element_text(size=20),
                                             text=element_text(size=20))
wt_ntx_ven_gg_lab
dev.copy(jpeg,"FC_nTG_nTX_ven.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX venules RFC = ",round(wt_ntx_ven_flow_change,digits=3)," +/- ",round(wt_ntx_ven_flow_change_ci,digits=3)," R2 = ",signif(summary(wt_ntx_ven_linfit)$adj.r.squared, 3)))
wt_ntx_ven_flow_change_tab = paste(round(wt_ntx_ven_flow_change,digits=3)," +/- ",round(wt_ntx_ven_flow_change_ci,digits=3))
# venule_regression tg ntx
airdata_tg_ntx_ven_mtchd<-subset(airdata_tg_ntx_mtchd,type=="v")
cdata_tg_ntx_ven_mtchd<-subset(cdata_tg_ntx_mtchd,type=="v")
tg_ntx_ven_ac <-data.frame(airdata_tg_ntx_ven_mtchd$pt,cdata_tg_ntx_ven_mtchd$pt)
tg_ntx_ven_linfit <- lm(cdata_tg_ntx_ven_mtchd.pt~airdata_tg_ntx_ven_mtchd.pt,data=tg_ntx_ven_ac)
tg_ntx_ven_ci <- confint(tg_ntx_ven_linfit, 'airdata_tg_ntx_ven_mtchd.pt', level=0.95)
tg_ntx_ven_ci_int <- (tg_ntx_ven_ci[2]-tg_ntx_ven_ci[1])/2
tg_ntx_ven_gg <- ggplotRegression(tg_ntx_ven_linfit)
tg_ntx_ven_linfit_sm <-summary(tg_ntx_ven_linfit)
tg_ntx_ven_linfit_std <- tg_ntx_ven_linfit_sm$coefficients[[4]]
tg_ntx_ven_flow_change <- 1/tg_ntx_ven_linfit$coefficients[2]
tg_ntx_ven_flow_change_ci <- (tg_ntx_ven_ci_int*tg_ntx_ven_flow_change)^2
# tg_ntx_ven_flow_change_std <- sqrt((tg_ntx_ven_linfit_std*tg_ntx_ven_flow_change/tg_ntx_ven_linfit$coefficients[2])^2)
tg_ntx_ven_gg_lab <- tg_ntx_ven_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("Tg - Venules")
tg_ntx_ven_gg_lab <- tg_ntx_ven_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
tg_ntx_ven_gg_lab <- tg_ntx_ven_gg_lab+theme(axis.text=element_text(size=20),
                                             text=element_text(size=20))
tg_ntx_ven_gg_lab
dev.copy(jpeg,"FC_TG_nTX_ven.jpeg")
dev.off()
#+ggtitle(paste("nTG nTX venules RFC = ",round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_change_ci,digits=3)," R2 = ",signif(summary(tg_ntx_ven_linfit)$adj.r.squared, 3)))
tg_ntx_ven_flow_change_tab = paste(round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_change_ci,digits=3))
# z for hypothesis testing of two linear regression
# Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_ven_flow_change = (tg_ntx_ven_linfit$coefficients[2]-wt_ntx_ven_linfit$coefficients[2])/
  (sqrt(tg_ntx_ven_linfit_std^2+wt_ntx_ven_linfit_std^2))
p_ntx_ven_flow_change=2*pnorm(-abs(z_ntx_ven_flow_change))
# /Ordinary regression (error in y only) - end
# -----------------------------------

## not fixed for ci p value nTG TX
air_not_c <- airdata_wt_tx[airdata_wt_tx$type!="c",]
c_not_c <- cdata_wt_tx[cdata_wt_tx$type!="c",]
data_c_air <- subset(airdata_wt_tx, type=="c")
data_c_co2 <- subset(cdata_wt_tx,  type=="c")
keys <- c("vessel_number", "depth")
data_c_air_match <- merge(data_c_air, data_c_co2,by=keys)
ves <- data_c_air_match$vessel_number
data_c_air_mtchd <- subset(data_c_air, vessel_number %in% ves)
data_c_co2_mtchd <- subset(data_c_co2, vessel_number %in% ves)
data_c_mtchd <- rbind(data_c_air_mtchd,data_c_co2_mtchd)
data.frame(table(data_c_mtchd$state))   # number of air and co2 should be the same
airdata_wt_tx_mtchd <- rbind(air_not_c,data_c_air_mtchd)
cdata_wt_tx_mtchd <- rbind(c_not_c,data_c_co2_mtchd)

# all_regression
wt_tx_ac <-data.frame(airdata_wt_tx_mtchd$pt,cdata_wt_tx_mtchd$pt)
wt_tx_linfit <- lm(cdata_wt_tx_mtchd.pt~airdata_wt_tx_mtchd.pt,data=wt_tx_ac)
wt_tx_gg <- ggplotRegression(wt_tx_linfit)
wt_tx_linfit_sm <-summary(wt_tx_linfit)
wt_tx_linfit_std <- wt_tx_linfit_sm$coefficients[[4]]
wt_tx_flow_change <- 1/wt_tx_linfit$coefficients[2]
wt_tx_flow_change_std <- sqrt((wt_tx_linfit_std*wt_tx_flow_change/wt_tx_linfit$coefficients[2])^2)
wt_tx_gg_lab <- wt_tx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle(paste("nTG TX rats RFC = ",round(wt_tx_flow_change,digits=3)," +/- ",round(wt_tx_flow_change_std,digits=3)," R2 = ",signif(summary(wt_tx_linfit)$adj.r.squared, 3)))
wt_tx_gg_lab <- wt_tx_gg_lab+xlim(4,9)+ylim(4,9)
wt_tx_gg_lab
dev.copy(jpeg,"FC_nTG_TX_rats.jpeg")
dev.off()

## TG TX
air_not_c <- airdata_tg_tx[airdata_tg_tx$type!="c",]
c_not_c <- cdata_tg_tx[cdata_tg_tx$type!="c",]
data_c_air <- subset(airdata_tg_tx, type=="c")
data_c_co2 <- subset(cdata_tg_tx,  type=="c")
keys <- c("vessel_number", "depth")
data_c_air_match <- merge(data_c_air, data_c_co2,by=keys)
ves <- data_c_air_match$vessel_number
data_c_air_mtchd <- subset(data_c_air, vessel_number %in% ves)
data_c_co2_mtchd <- subset(data_c_co2, vessel_number %in% ves)
data_c_mtchd <- rbind(data_c_air_mtchd,data_c_co2_mtchd)
data.frame(table(data_c_mtchd$state))   # number of air and co2 should be the same
airdata_tg_tx_mtchd <- rbind(air_not_c,data_c_air_mtchd)
cdata_tg_tx_mtchd <- rbind(c_not_c,data_c_co2_mtchd)

# all_regression
tg_tx_ac <-data.frame(airdata_tg_tx_mtchd$pt,cdata_tg_tx_mtchd$pt)
tg_tx_linfit <- lm(cdata_tg_tx_mtchd.pt~airdata_tg_tx_mtchd.pt,data=tg_tx_ac)
tg_tx_gg <- ggplotRegression(tg_tx_linfit)
tg_tx_linfit_sm <-summary(tg_tx_linfit)
tg_tx_linfit_std <- tg_tx_linfit_sm$coefficients[[4]]
tg_tx_flow_change <- 1/tg_tx_linfit$coefficients[2]
tg_tx_flow_change_std <- sqrt((tg_tx_linfit_std*tg_tx_flow_change/tg_tx_linfit$coefficients[2])^2)
tg_tx_gg_lab <- tg_tx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle(paste("TG TX rats RFC = ",round(tg_tx_flow_change,digits=3)," +/- ",round(tg_tx_flow_change_std,digits=3)," R2 = ",signif(summary(tg_tx_linfit)$adj.r.squared, 3)))
tg_tx_gg_lab <- tg_tx_gg_lab+xlim(6.5,13)+ylim(6.5,13)
tg_tx_gg_lab
dev.copy(jpeg,"FC_TG_TX_rats.jpeg")
dev.off()


#TEMPLATE # Ordinary least square (minimizing difference in Y only) linear fitting
# wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
# wt_ntx_linfit <- lm(cdata_wt_ntx_mtchd.pt~airdata_wt_ntx_mtchd.pt,data=wt_ntx_ac)
# wt_ntx_ci <- confint(wt_ntx_linfit, 'airdata_wt_ntx_mtchd.pt', level=0.95)
# wt_ntx_ci_int <- (wt_ntx_ci[2]-wt_ntx_ci[1])/2
# wt_ntx_gg <- ggplotRegression(wt_ntx_linfit)
# wt_ntx_linfit_sm <-summary(wt_ntx_linfit)
# wt_ntx_linfit_std <- wt_ntx_linfit_sm$coefficients[[4]]
# wt_ntx_flow_change <- 1/wt_ntx_linfit$coefficients[2]
# wt_ntx_flow_change_ci <- (wt_ntx_ci_int*wt_ntx_flow_change)^2
# # wt_ntx_flow_change_std <- sqrt((wt_ntx_linfit_std*wt_ntx_flow_change/wt_ntx_linfit$coefficients[2])^2)
# wt_ntx_gg_lab <- wt_ntx_gg+labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg")
# wt_ntx_gg_lab <- wt_ntx_gg_lab+xlim(4,12)+ylim(4,12)+geom_abline(slope=1,intercept=0,colour="blue")
# wt_ntx_gg_lab <- wt_ntx_gg_lab+theme(axis.text=element_text(size=20),
#                                      text=element_text(size=20))
# wt_ntx_gg_lab
#TEMPLATE


all_ <- lm(cdata.pt~airdata.pt,data=aa)
ggplotRegression(linfit)

flow_change <- 1/linfit$coefficients[2]

means<-aggregate(PT~State+Vessel_Number, tr1, mean)

# pct100_60sec_regression
pct100_60sec_linfit <- lm(pct100_60sec_cdata.pt~pct100_60sec_airdata.pt,data=pct100_60sec_ac)
pct100_60sec_gg <- ggplotRegression(pct100_60sec_linfit)
pct100_60sec_linfit_sm <-summary(pct100_60sec_linfit)
pct100_60sec_linfit_std <- pct100_60sec_linfit_sm$coefficients[[4]]
pct100_60sec_flow_change <- 1/pct100_60sec_linfit$coefficients[2]
pct100_60sec_flow_change_std <- sqrt((pct100_60sec_linfit_std*pct100_60sec_flow_change/pct100_60sec_linfit$coefficients[2])^2)
pct100_60sec_gg_lab <- pct100_60sec_gg+labs(y="TTP during CO2 breating [s]")+labs(x="TTP during air breating [s]")+ggtitle(paste("ALL - 10% challenge (60sec) RFC = ",round(pct100_60sec_flow_change,digits=3)," +/- ",round(pct100_60sec_flow_change_std,digits=3)))
dev.copy(jpeg,"FC_pct100_60sec_all.jpeg")
dev.off()



### box plot showing the variance of TTP in each state from each model 
data_46 = subset(data,id==46)
data_48 = subset(data,id==48)
data_40 = subset(data,id==40)
data_51 = subset(data,id==51)

airdata=subset(data,data$state=="air")
qplot(airdata$vessel_number,airdata$pt,colour = airdata$depth, shape = airdata$)
cdata=subset(data,data$state=="co2")
qplot(cdata$id,cdata$pt,colour = cdata$dept) 

qplot(data_46$vessel_number,data_46$pt,colour = data_46$depth,shape = data_46$state)
qplot(data_48$vessel_number,data_48$pt,colour = data_48$depth,shape = data_48$state)
qplot(data_40$vessel_number,data_40$pt,colour = data_40$depth,shape = data_40$state)
qplot(data_51$vessel_number,data_51$pt,colour = data_51$depth,shape = data_51$state)
data_48_depth_air = (mean(data[data_48$state=="air"&data_48$depth==150,]$pt) -
              mean(data[data_48$state=="air"&data_48$depth==200,]$pt))/
              mean(data[data_48$state=="air"&data_48$depth==150,]$pt)
data_48_depth_air
data_48_depth_co2 = (mean(data[data_48$state=="co2"&data_48$depth==150,]$pt) -
                   mean(data[data_48$state=="co2"&data_48$depth==200,]$pt))/
                   mean(data[data_48$state=="co2"&data_48$depth==150,]$pt)
data_48_depth_co2
data_48_150_airtoco2 = (mean(data[data_48$state=="co2"&data_48$depth==150,]$pt) -
                        mean(data[data_48$state=="air"&data_48$depth==150,]$pt))/
                        mean(data[data_48$state=="air"&data_48$depth==150,]$pt)
data_48_150_airtoco2
data_48_200_airtoco2 = (mean(data[data_48$state=="co2"&data_48$depth==200,]$pt) -
                          mean(data[data_48$state=="air"&data_48$depth==200,]$pt))/
  mean(data[data_48$state=="air"&data_48$depth==200,]$pt)
data_48_200_airtoco2
data_48_airtoco2 = (mean(data[data_48$state=="co2",]$pt) -
                          mean(data[data_48$state=="air",]$pt))/
  mean(data[data_48$state=="air",]$pt)
data_48_airtoco2

qplot(data_48$vessel_number,data_48$pt,colour = data_48$depth,shape = data_48$state)
cdata_46=subset(data_46,data_46$state=="co2")
qplot(cdata_46$vessel_number,cdata_46$pt,colour = cdata_46$id) 

data_summary <- summarySE(data, measurevar="pt", groupvars=c("gntrt","state"))
data_summ <- data_summary
data_summ
data_summ$gntrt = factor(data_summ$gntrt)
# TTP box
ggplot(data_summ, aes(x=gntrt, y=pt, fill=state)) +
  geom_bar(position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(ymin=pt, ymax=pt+sd),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Genotype & Treatment") +
  ylab("Time to peak (s)") +
  scale_fill_hue(name="State")+
  ggtitle("The effect of co2 breathing on the time to peak")+
  scale_x_discrete("Genotype & Treatment",labels=c("tg&ntx" = "Tg",
                                                   "tg&tx" = "Tg LNAME",
                                                   "wt&ntx" = "nTg",
                                                   "wt&tx" = "nTg LNAME"))
dev.copy(jpeg,"TTP_all_tgVStx.jpeg")
dev.off()
cov_tg_ntx_air<-data_summ[data_summ$gntrt=="tg&ntx"&data_summ$state=="air",]$sd/
  data_summ[data_summ$gntrt=="tg&ntx"&data_summ$state=="air",]$pt
cov_tg_ntx_air
cov_tg_ntx_co2<-data_summ[data_summ$gntrt=="tg&ntx"&data_summ$state=="co2",]$sd/
  data_summ[data_summ$gntrt=="tg&ntx"&data_summ$state=="co2",]$pt
cov_tg_ntx_co2
cov_tg_tx_air<-data_summ[data_summ$gntrt=="tg&tx"&data_summ$state=="air",]$sd/
  data_summ[data_summ$gntrt=="tg&tx"&data_summ$state=="air",]$pt
cov_tg_tx_air
cov_tg_tx_co2<-data_summ[data_summ$gntrt=="tg&tx"&data_summ$state=="co2",]$sd/
  data_summ[data_summ$gntrt=="tg&tx"&data_summ$state=="co2",]$pt
cov_tg_tx_co2
cov_wt_ntx_air<-data_summ[data_summ$gntrt=="wt&ntx"&data_summ$state=="air",]$sd/
  data_summ[data_summ$gntrt=="wt&ntx"&data_summ$state=="air",]$pt
cov_wt_ntx_air
cov_wt_ntx_co2<-data_summ[data_summ$gntrt=="wt&ntx"&data_summ$state=="co2",]$sd/
  data_summ[data_summ$gntrt=="wt&ntx"&data_summ$state=="co2",]$pt
cov_wt_ntx_co2
cov_wt_tx_air<-data_summ[data_summ$gntrt=="wt&tx"&data_summ$state=="air",]$sd/
  data_summ[data_summ$gntrt=="wt&tx"&data_summ$state=="air",]$pt
cov_wt_tx_air
cov_wt_tx_co2<-data_summ[data_summ$gntrt=="wt&tx"&data_summ$state=="co2",]$sd/
  data_summ[data_summ$gntrt=="wt&tx"&data_summ$state=="co2",]$pt
cov_wt_tx_co2
cov_tg_ntx<-c(cov_tg_ntx_air,cov_tg_ntx_co2)
cov_tg_tx<-c(cov_tg_tx_air,cov_tg_tx_co2)
cov_wt_ntx<-c(cov_wt_ntx_air,cov_wt_ntx_co2)
cov_wt_tx<-c(cov_wt_tx_air,cov_wt_tx_co2)

cov <- data.frame(cov_tg_ntx,cov_tg_tx,cov_wt_ntx,cov_wt_tx)
colnames(cov)<-c("Tg","Tg LNAME","nTg","nTg LNAME")
rownames(cov)<-c("air","co2")
cov
data_150 = subset(data,data$depth==150)
data_summary_150 <- summarySE(data_150, measurevar="pt", groupvars=c("gntrt","state"))
data_summ_150 <- data_summary_150
data_summ_150
# TTP box
ggplot(data_summ_150, aes(x=gntrt, y=pt, fill=state)) +
  geom_bar(position=position_dodge(), stat="identity",
           color="black",
           size=.3) +
  geom_errorbar(aes(ymin=pt, ymax=pt+sd),
                size=.3,
                width=.2,
                position=position_dodge(.9)) +
  xlab("Genotype & Treatment") +
  ylab("Time to peak (s)") +
  scale_fill_hue(name="State")+
  ggtitle("The effect of co2 breathing on the time to peak, depth=150um")+
  scale_x_discrete("Genotype & Treatment",labels=c("tg&ntx" = "Tg",
                                                   "tg&tx" = "Tg LNAME",
                                                   "wt&ntx" = "nTg",
                                                   "wt&tx" = "nTg LNAME"))
dev.copy(jpeg,"TTP_150_tgVStx.jpeg")
dev.off()
cov_tg_ntx_air<-data_summ_150[data_summ_150$gntrt=="tg&ntx"&data_summ_150$state=="air",]$sd/
  data_summ_150[data_summ_150$gntrt=="tg&ntx"&data_summ_150$state=="air",]$pt
cov_tg_ntx_air
cov_tg_ntx_co2<-data_summ_150[data_summ_150$gntrt=="tg&ntx"&data_summ_150$state=="co2",]$sd/
  data_summ_150[data_summ_150$gntrt=="tg&ntx"&data_summ_150$state=="co2",]$pt
cov_tg_ntx_co2
cov_tg_tx_air<-data_summ_150[data_summ_150$gntrt=="tg&tx"&data_summ_150$state=="air",]$sd/
  data_summ_150[data_summ_150$gntrt=="tg&tx"&data_summ_150$state=="air",]$pt
cov_tg_tx_air
cov_tg_tx_co2<-data_summ_150[data_summ_150$gntrt=="tg&tx"&data_summ_150$state=="co2",]$sd/
  data_summ_150[data_summ_150$gntrt=="tg&tx"&data_summ_150$state=="co2",]$pt
cov_tg_tx_co2
cov_wt_ntx_air<-data_summ_150[data_summ_150$gntrt=="wt&ntx"&data_summ_150$state=="air",]$sd/
  data_summ_150[data_summ_150$gntrt=="wt&ntx"&data_summ_150$state=="air",]$pt
cov_wt_ntx_air
cov_wt_ntx_co2<-data_summ_150[data_summ_150$gntrt=="wt&ntx"&data_summ_150$state=="co2",]$sd/
  data_summ_150[data_summ_150$gntrt=="wt&ntx"&data_summ_150$state=="co2",]$pt
cov_wt_ntx_co2
cov_wt_tx_air<-data_summ_150[data_summ_150$gntrt=="wt&tx"&data_summ_150$state=="air",]$sd/
  data_summ_150[data_summ_150$gntrt=="wt&tx"&data_summ_150$state=="air",]$pt
cov_wt_tx_air
cov_wt_tx_co2<-data_summ_150[data_summ_150$gntrt=="wt&tx"&data_summ_150$state=="co2",]$sd/
  data_summ_150[data_summ_150$gntrt=="wt&tx"&data_summ_150$state=="co2",]$pt
cov_wt_tx_co2
cov_tg_ntx<-c(cov_tg_ntx_air,cov_tg_ntx_co2)
cov_tg_tx<-c(cov_tg_tx_air,cov_tg_tx_co2)
cov_wt_ntx<-c(cov_wt_ntx_air,cov_wt_ntx_co2)
cov_wt_tx<-c(cov_wt_tx_air,cov_wt_tx_co2)

cov <- data.frame(cov_tg_ntx,cov_tg_tx,cov_wt_ntx,cov_wt_tx)
colnames(cov)<-c("Tg","Tg LNAME","nTg","nTg LNAME")
rownames(cov)<-c("air","co2")
cov

# Slope
# lme MODEL: flow ~ state + gn + (1|subject) + (1|fl accumation) + e
