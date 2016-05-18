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
data<-as.data.frame(read.csv("bolus_AD_Rat_all_OLrej_rev5_overlap.csv",header=T))
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

airdata_all<-subset(data_matched,data_matched$state=="air")
cdata_all<-subset(data_matched,data_matched$state=="co2")
delta_all_pt <- airdata_all$pt - cdata_all$pt
delta_all_p_pt <- (airdata_all$pt - cdata_all$pt)/airdata_all$pt*100
airdata_all$delta_pt <- delta_all_pt
airdata_all$delta_p_pt <- delta_all_p_pt
airdata_all <- airdata_all[order(rev(airdata_all$gntrt),airdata_all$id,airdata_all$vessel_number),] #
airdata_all <- subset(airdata_all,id!=74) # excluded due to imaging depth shift
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
data_summary <- summarySE(airdata_all_all, measurevar="delta_p_pt", groupvars=c("type","gntrt"))
data_summ <- data_summary
data_summ
# 4 bars in each group
desired_order <- c("wt&ntx","wt&tx","tg&ntx","tg&tx")
data_summ$gntrt <- factor(as.character(data_summ$gntrt), levels=desired_order)
data_summ <- data_summ[order(data_summ$type,data_summ$gntrt),]
temp <- data_summ[5:8,]
data_summ[5:8,] <- data_summ[13:16,]
data_summ[13:16,] <- temp
x.seq <- c(1,2,3,4,6.5,7.5,8.5,9.5,11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5)
data_summ$x <- x.seq
data_summ
aa<-ggplot(data_summ, aes(x=x, y=delta_p_pt)) +
  geom_bar(aes(fill=gntrt),position=position_dodge(), stat="identity",
           color="black",
           width=1,
           size=1) +
  geom_errorbar(aes(fill=gntrt,ymin=delta_p_pt, ymax=delta_p_pt+se),
                size=.7,
                width=.4,
                position=position_dodge(.9)) +
  xlab("") +  ylab("Vascular Reactivity [%]") + ylim(0,25)+
  scale_fill_manual(name="Genotype&Treatment",
                    values=c("darkolivegreen2","olivedrab4","darkorange","olivedrab3"),labels=c("nTg","nTg+LNAME","Tg","Tg+LNAME"))+
  scale_x_continuous(breaks=c(sum(x.seq[1:4])/4, sum(x.seq[5:8])/4, sum(x.seq[9:12])/4, sum(x.seq[13:16])/4),
                     label=c("All Vessels","Arterioles","Capillaries","Venules"))+theme_gray()+
  theme(axis.text=element_text(size=30),
        axis.title.y=element_text(vjust=5),
        text=element_text(size=30),
        legend.position="bottom")
aa
label1.df <- data.frame(x = c(x.seq[3],x.seq[2],x.seq[7],x.seq[6],x.seq[15],x.seq[14])-0.1,
                        delta_p_pt = c(19.5,18.5,20.5,19.5,22.5,21.5))
data1 <- data.frame(x = c(1,3), y = c(19.5,19.5),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(1,2), y = c(18.5,18.5),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(6.5,8.5), y = c(20.5,20.5),
                    type=NA, delta_pt=NA)
data4 <- data.frame(x = c(6.5,7.5), y = c(19.5,19.5),
                    type=NA, delta_pt=NA)
data5 <- data.frame(x = c(16.5,18.5), y = c(22.5,22.5),
                    type=NA, delta_pt=NA)
data6 <- data.frame(x = c(16.5,17.5), y = c(21.5,21.5),
                    type=NA, delta_pt=NA)
bb<-aa + geom_path(data = data1, aes(x = x, y = y),size=1)+
  geom_path(data = data2, aes(x = x, y = y),size=1)+
  geom_path(data = data3, aes(x = x, y = y),size=1)+
  geom_path(data = data4, aes(x = x, y = y),size=1)+
  geom_path(data = data5, aes(x = x, y = y),size=1)+
  geom_path(data = data6, aes(x = x, y = y),size=1)+
  geom_text(data=label1.df,label="*",size=13)
bb
dev.copy(jpeg,"LNAME/vascular_reac_LNAME.jpeg", width=1000, height=800)
dev.off()
airdata_all_all_all <- subset(airdata_all_all,type=="All")
wilcox.test(airdata_all_all_all[airdata_all_all_all$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_all[airdata_all_all_all$gntrt=="wt&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_all[airdata_all_all_all$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_all[airdata_all_all_all$gntrt=="tg&ntx",]$delta_p_pt)
wilcox.test(airdata_all_all_all[airdata_all_all_all$gntrt=="wt&tx",]$delta_p_pt,airdata_all_all_all[airdata_all_all_all$gntrt=="tg&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_all[airdata_all_all_all$gntrt=="tg&ntx",]$delta_p_pt,airdata_all_all_all[airdata_all_all_all$gntrt=="tg&tx",]$delta_p_pt)
airdata_all_all_art <- subset(airdata_all_all,type=="a")
wilcox.test(airdata_all_all_art[airdata_all_all_art$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_art[airdata_all_all_art$gntrt=="wt&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_art[airdata_all_all_art$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_art[airdata_all_all_art$gntrt=="tg&ntx",]$delta_p_pt)
wilcox.test(airdata_all_all_art[airdata_all_all_art$gntrt=="wt&tx",]$delta_p_pt,airdata_all_all_art[airdata_all_all_art$gntrt=="tg&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_art[airdata_all_all_art$gntrt=="tg&ntx",]$delta_p_pt,airdata_all_all_art[airdata_all_all_art$gntrt=="tg&tx",]$delta_p_pt)
airdata_all_all_cap <- subset(airdata_all_all,type=="c")
wilcox.test(airdata_all_all_cap[airdata_all_all_cap$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_cap[airdata_all_all_cap$gntrt=="wt&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_cap[airdata_all_all_cap$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_cap[airdata_all_all_cap$gntrt=="tg&ntx",]$delta_p_pt)
wilcox.test(airdata_all_all_cap[airdata_all_all_cap$gntrt=="wt&tx",]$delta_p_pt,airdata_all_all_cap[airdata_all_all_cap$gntrt=="tg&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_cap[airdata_all_all_cap$gntrt=="tg&ntx",]$delta_p_pt,airdata_all_all_cap[airdata_all_all_cap$gntrt=="tg&tx",]$delta_p_pt)
airdata_all_all_ven <- subset(airdata_all_all,type=="v")
wilcox.test(airdata_all_all_ven[airdata_all_all_ven$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_ven[airdata_all_all_ven$gntrt=="wt&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_ven[airdata_all_all_ven$gntrt=="wt&ntx",]$delta_p_pt,airdata_all_all_ven[airdata_all_all_ven$gntrt=="tg&ntx",]$delta_p_pt)
wilcox.test(airdata_all_all_ven[airdata_all_all_ven$gntrt=="wt&tx",]$delta_p_pt,airdata_all_all_ven[airdata_all_all_ven$gntrt=="tg&tx",]$delta_p_pt)
wilcox.test(airdata_all_all_ven[airdata_all_all_ven$gntrt=="tg&ntx",]$delta_p_pt,airdata_all_all_ven[airdata_all_all_ven$gntrt=="tg&tx",]$delta_p_pt) 


# REGRESSIONS
# "cowplot" library is needed for regression plot with background shading but it messes up the default ggplot 
# theme...either add "theme_gray()" to ggplot syntax or restart R and plot without "cowplot" library
library(cowplot)
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
ggplotRegressionPCAcentered <- function (df,x_left,x_right) {
  require(ggplot2)
  ggplot(df,aes_string(names(df)[3],names(df)[4]))+geom_point()+
    geom_line(aes(x=TP_air_centered,y=cu),colour="darkgrey",alpha="0.5")+
    geom_line(aes(x=TP_air_centered,y=cl),colour="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cl,ymax=cu),fill="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cu,ymax=cl),fill="darkgrey",alpha="0.5")+
    stat_smooth(data=subset(df,TP_co2_centered>x_left&TP_co2_centered<x_right),
                method = f, se = TRUE, colour = 'red', size=1.5, formula=y~x)
}
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

# CI for Orthogonal Least Square fitting of a line: http://stackoverflow.com/questions/31057192/resampling-not-producing-expected-result-of-principal-component-analysis
## nTG nTX - matching air&co2
regPrepCentered <- function (df,genotrt) {
  air_co2 <-data.frame(df[df$gntrt==genotrt&df$state=="air",]$pt,df[df$gntrt==genotrt&df$state=="co2",]$pt)
  colnames(air_co2) <- c("air","co2")
  air_co2$TP_air_centered <- air_co2$air - mean(air_co2$air,na.rm=TRUE)
  air_co2$TP_co2_centered <- air_co2$co2 - mean(air_co2$co2,na.rm=TRUE)
  air_co2_pca <- deming(air_co2[,4]~air_co2[,3])
  air_co2_flow_change <- 1/air_co2_pca$coefficients[2]
  air_co2_ci_int <- (air_co2_pca$ci[2,2]-air_co2_pca$ci[2,1])/2
  air_co2_flow_change_ci <- air_co2_ci_int*(air_co2_flow_change)^2
  air_co2_slp_se <- sqrt(diag(air_co2_pca$variance))[2]                     # se for slope
  air_co2_flow_se <- air_co2_slp_se/air_co2_pca$coefficients[2]^2            # se for flow change
  air_co2_slp_int_se <- sqrt(diag(air_co2_pca$variance))[1]                 # se for original slope's intercept
  air_co2_flow_int_se <- air_co2_slp_int_se/air_co2_pca$coefficients[1]^2    # se for flow change
  air_co2$se_u <- (air_co2_pca$coefficients[2]+air_co2_slp_se)*air_co2$TP_air_centered+air_co2_pca$coeff[1]-air_co2_slp_int_se
  air_co2$se_l <- (air_co2_pca$coefficients[2]-air_co2_slp_se)*air_co2$TP_air_centered+air_co2_pca$coeff[1]+air_co2_slp_int_se
  air_co2$cu <- air_co2_pca$ci[2,2]*air_co2$TP_air_centered+air_co2_pca$ci[1,1]
  air_co2$cl <- air_co2_pca$ci[2,1]*air_co2$TP_air_centered+air_co2_pca$ci[1,2]
  return(air_co2)
}
regPrep <- function (df,genotrt) {
  air_co2 <-data.frame(df[df$gntrt==genotrt&df$state=="air",]$pt,df[df$gntrt==genotrt&df$state=="co2",]$pt)
  colnames(air_co2) <- c("air","co2")
  air_co2_pca <- deming(air_co2[,2]~air_co2[,1])
  air_co2_flow_change <- 1/air_co2_pca$coefficients[2]
  air_co2_ci_int <- (air_co2_pca$ci[2,2]-air_co2_pca$ci[2,1])/2              # ci interval for slope
  air_co2_flow_change_ci <- air_co2_ci_int*(air_co2_flow_change)^2           # ci interval for flow change
  air_co2_slp_se <- sqrt(diag(air_co2_pca$variance))[2]                     # se for slope
  air_co2_flow_se <- air_co2_slp_se/air_co2_pca$coefficients[2]^2            # se for flow change
  air_co2_slp_int_se <- sqrt(diag(air_co2_pca$variance))[1]                 # se for original slope's intercept
  air_co2_flow_int_se <- air_co2_slp_int_se/air_co2_pca$coefficients[1]^2    # se for flow change
  air_co2$se_u <- (air_co2_pca$coefficients[2]+air_co2_slp_se)*air_co2$air+air_co2_pca$coeff[1]-air_co2_slp_int_se
  air_co2$se_l <- (air_co2_pca$coefficients[2]-air_co2_slp_se)*air_co2$air+air_co2_pca$coeff[1]+air_co2_slp_int_se
  air_co2$cu <- air_co2_pca$ci[2,2]*air_co2$air+air_co2_pca$ci[1,1]
  air_co2$cl <- air_co2_pca$ci[2,1]*air_co2$air+air_co2_pca$ci[1,2]
  return(air_co2)
}
flowChange <- function (df,genotrt) {
  air_co2 <-data.frame(df[df$gntrt==genotrt&df$state=="air",]$pt,df[df$gntrt==genotrt&df$state=="co2",]$pt)
  colnames(air_co2) <- c("air","co2")
  air_co2_pca <- deming(air_co2[,2]~air_co2[,1])
  air_co2_flow_change <- 1/air_co2_pca$coefficients[2]
  air_co2_ci_int <- (air_co2_pca$ci[2,2]-air_co2_pca$ci[2,1])/2              # ci interval for slope
  air_co2_flow_change_ci <- air_co2_ci_int*(air_co2_flow_change)^2           # ci interval for flow change
  air_co2_slp_se <- sqrt(diag(air_co2_pca$variance))[2]                     # se for slope
  air_co2_flow_se <- air_co2_slp_se/air_co2_pca$coefficients[2]^2            # se for flow change
  air_co2_flow_change_tab = paste(round(air_co2_flow_change,digits=3)," +/- ",round(air_co2_flow_se,digits=3))
  return(air_co2_flow_change_tab)
}
sigP <- function (df,genotrt1,genotrt2) {
  air_co2_1 <-data.frame(df[df$gntrt==genotrt1&df$state=="air",]$pt,df[df$gntrt==genotrt1&df$state=="co2",]$pt)
  colnames(air_co2_1) <- c("air","co2")
  air_co2_1_pca <- deming(air_co2_1[,2]~air_co2_1[,1])
  air_co2_1_flow_change <- 1/air_co2_1_pca$coefficients[2]
  air_co2_1_ci_int <- (air_co2_1_pca$ci[2,2]-air_co2_1_pca$ci[2,1])/2              # ci interval for slope
  air_co2_1_flow_change_ci <- air_co2_1_ci_int*(air_co2_1_flow_change)^2           # ci interval for flow change
  air_co2_1_slp_se <- sqrt(diag(air_co2_1_pca$variance))[2]                     # se for slope
  air_co2_1_flow_se <- air_co2_1_slp_se/air_co2_1_pca$coefficients[2]^2            # se for flow change
  air_co2_1_slp_int_se <- sqrt(diag(air_co2_1_pca$variance))[1]                 # se for original slope's intercept
  air_co2_1_flow_int_se <- air_co2_1_slp_int_se/air_co2_1_pca$coefficients[1]^2    # se for flow change
  
  air_co2_2 <-data.frame(df[df$gntrt==genotrt2&df$state=="air",]$pt,df[df$gntrt==genotrt2&df$state=="co2",]$pt)
  colnames(air_co2_2) <- c("air","co2")
  air_co2_2_pca <- deming(air_co2_2[,2]~air_co2_2[,1])
  air_co2_2_flow_change <- 1/air_co2_2_pca$coefficients[2]
  air_co2_2_ci_int <- (air_co2_2_pca$ci[2,2]-air_co2_2_pca$ci[2,1])/2              # ci interval for slope
  air_co2_2_flow_change_ci <- air_co2_2_ci_int*(air_co2_2_flow_change)^2           # ci interval for flow change
  air_co2_2_slp_se <- sqrt(diag(air_co2_2_pca$variance))[2]                     # se for slope
  air_co2_2_flow_se <- air_co2_2_slp_se/air_co2_2_pca$coefficients[2]^2            # se for flow change
  air_co2_2_slp_int_se <- sqrt(diag(air_co2_2_pca$variance))[1]                 # se for original slope's intercept
  air_co2_2_flow_int_se <- air_co2_2_slp_int_se/air_co2_2_pca$coefficients[1]^2    # se for flow change
  
  z_air_co2_flow_change = (air_co2_1_pca$coefficient[2]-air_co2_2_pca$coefficient[2])/
    (sqrt(air_co2_1_slp_se^2+air_co2_2_slp_se^2))
  p_air_co2_flow_change=2*pnorm(-abs(z_air_co2_flow_change))
  return(p_air_co2_flow_change)
}
# ORTHOGONAL REGRESSION - CENTERED - ALL VESSELS
wt_ntx_ac <- regPrepCentered(data_matched,"wt&ntx")
wt_tx_ac <- regPrepCentered(data_matched,"wt&tx")
tg_ntx_ac <- regPrepCentered(data_matched,"tg&ntx")
tg_tx_ac <- regPrepCentered(data_matched,"tg&tx")

# all_regression wt ntx     CENTRED ALL
wt_ntx_plot <- ggplotRegressionPCAcentered(wt_ntx_ac,-3,3.5)
wt_ntx_plot_cnt <- wt_ntx_plot+xlim(-3,3.5)+ylim(-3,3.5)+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("nTg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
wt_ntx_plot_cnt
wt_ntx_flow_change_tab = flowChange(data_matched,"wt&ntx")
# all_regression wt tx     CENTRED ALL
wt_tx_plot <- ggplotRegressionPCAcentered(wt_tx_ac,-3,4.5)
wt_tx_plot_cnt <- wt_tx_plot+xlim(-3,4.5)+ylim(-3,4.5)+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("nTg+LNAME - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
wt_tx_plot_cnt
wt_tx_flow_change_tab = flowChange(data_matched,"wt&tx")
# all_regression tg ntx     CENTRED ALL
tg_ntx_plot <- ggplotRegressionPCAcentered(tg_ntx_ac,-2.5,4.5)
tg_ntx_plot_cnt <- tg_ntx_plot+xlim(-2.5,4.5)+ylim(-2.5,4.5)+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("Tg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
tg_ntx_plot_cnt
tg_ntx_flow_change_tab = flowChange(data_matched,"tg&ntx")
# all_regression tg tx     CENTRED ALL
tg_tx_plot <- ggplotRegressionPCAcentered(tg_tx_ac,-3,5.5)
tg_tx_plot_cnt <- tg_tx_plot+xlim(-3,5.5)+ylim(-3,5.5)+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("Tg+LNAME - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
tg_tx_plot_cnt
tg_tx_flow_change_tab = flowChange(data_matched,"tg&tx")

# Combined plot of all vessels of nTg and all vessels of Tg     CENTERED
ggplotRegressionPCAcombined <- function (df1,df2,df3,df4,colorVar1,colorVar2,colorVar3,colorVar4,
                                         x_left1,x_right1,x_left2,x_right2,x_left3,x_right3,x_left4,x_right4) {
  require(ggplot2)
  # create a data frame with the total number of rows equal to all datapoint (tg+ntg+tgNAME+ntgLNAME) 
  new_df <- data.frame(gntrt=c(rep("Tg",nrow(df1)), rep("nTg",nrow(df2)),rep("Tg+LNAME",nrow(df3)), rep("nTg+LNAME",nrow(df4))))
  new_df <- cbind(new_df,rbind(df1,df2,df3,df4)) # bind the input data as new columns
  new_df <- cbind(new_df,data.frame(matrix(ncol=6,nrow=length(new_df[,1])))) # bind 6 more columns for division of flow increase and flow decrease regions
  left_lim <- floor(min(new_df$TP_air_centered,new_df$TP_co2_centered)/0.5)*0.5
  right_lim <- ceiling(max(new_df$TP_air_centered,new_df$TP_co2_centered)/0.5)*0.5
  lim <- max(abs(left_lim),right_lim)
  colnames(new_df)[(ncol(new_df)-5):ncol(new_df)] <- c("x","ydiag","ydiagneg","yhoriz","yvert1","yvert2")
  new_df[,(ncol(new_df)-5)] <- seq(-6,6,length=length(new_df[,1]))
  new_df[,(ncol(new_df)-4)] <- seq(-6,6,length=length(new_df[,1]))
  new_df[,(ncol(new_df)-3)] <- seq(6,-6,length=length(new_df[,1]))
  new_df[,(ncol(new_df)-2)] <- rep(0,length(new_df[,1]))
  # below: couldn't figure out how to creat y values for vertical line while following x values, so just created 
  # infinitely small and large values.
  new_df[,(ncol(new_df)-1)] <- c(rep(-9999999,floor(length(new_df[,1])/2)),rep(9999999,ceiling(length(new_df[,1])/2)))
  new_df[,ncol(new_df)] <- c(rep(9999999,floor(length(new_df[,1])/2)),rep(-9999999,ceiling(length(new_df[,1])/2)))
  
  ggplot(new_df,aes_string(names(new_df)[4],names(new_df)[5]))+geom_point(aes(color=gntrt),size=2,alpha=0.4)+
    geom_line(data=new_df,aes(x=x,y=ydiag),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=new_df,aes(x=x,y=ydiagneg),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=new_df,aes(x=x,y=yhoriz),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=new_df,aes(x=x,y=yvert1),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=new_df,aes(x=TP_air_centered,y=cu,group=gntrt),colour="darkgrey",alpha=0.5)+
    geom_line(data=new_df,aes(x=TP_air_centered,y=cl,group=gntrt),colour="darkgrey",alpha=0.5)+
    geom_ribbon(data=new_df,aes(x=x,ymin=ydiagneg,ymax=ydiag),fill="olivedrab3",alpha=0.3)+
    geom_ribbon(data=new_df,aes(x=x,ymin=ydiag,ymax=yvert1),fill="darkorange2",alpha=0.3)+
    geom_ribbon(data=new_df,aes(x=x,ymin=ydiagneg,ymax=yvert2),fill="darkorange2",alpha=0.3)+
    geom_ribbon(data=new_df,aes(ymin=cl,ymax=cu,group=gntrt),fill="darkgrey",alpha=0.6)+
    scale_color_manual(name="Genotype&Treatment",values=c(colorVar4,colorVar3,colorVar2,colorVar1),
                       labels=c("nTg+LNAME","Tg+LNAME","nTg","Tg"))+
    #     geom_point(data=df2,color=colorVar2,size=2,alpha=0.4)+
    #     geom_line(data=df2,aes(x=air,y=cu),colour="darkgrey",alpha="0.5")+
    #     geom_line(data=df2,aes(x=air,y=cl),colour="darkgrey",alpha="0.5")+
    #     geom_ribbon(data=df2,aes(ymin=cl,ymax=cu),fill="darkgrey",alpha=0.4)+
    #     geom_ribbon(data=df2,aes(ymin=cu,ymax=cl),fill="darkgrey",alpha=0.4)+
    stat_smooth(data=subset(df1,TP_air_centered>x_left1&TP_air_centered<x_right1),
                method = f, se = TRUE, colour = colorVar1, size=1.5, formula=y~x)+
    stat_smooth(data=subset(df2,TP_air_centered>x_left2&TP_air_centered<x_right2),
                method = f, se = TRUE, colour = colorVar2, size=1.5, formula=y~x)+
    stat_smooth(data=subset(df3,TP_air_centered>x_left3&TP_air_centered<x_right3),
                method = f, se = TRUE, colour = colorVar3, size=1.5, formula=y~x)+
    stat_smooth(data=subset(df4,TP_air_centered>x_left4&TP_air_centered<x_right4),
                method = f, se = TRUE, colour = colorVar4, size=1.5, formula=y~x)+
    coord_cartesian(xlim=c(-lim,lim),ylim=c(-lim,lim))
  
  #     geom_ribbon(data=new_df,aes(ymin=yhoriz,ymax=ydiag),fill="olivedrab4",alpha=0.2)+
  #     geom_ribbon(data=subset(new_df,x<0&x>-3),aes(ymin=ydiag,ymax=ydiagneg),fill="olivedrab4",alpha=0.2)+
}
combined_plot <- ggplotRegressionPCAcombined(tg_ntx_ac,wt_ntx_ac,tg_tx_ac,wt_tx_ac,
                                             "darkorange","olivedrab2","darkorange3","olivedrab4",
                                             -2.5,4.5,-3,3.5,-3,5.5,-2.5,4.5)
combined_plot

comb <- combined_plot+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+theme_gray()+
  theme(axis.text=element_text(size=30),text=element_text(size=30),
        legend.position="none")
comb
dev.copy(jpeg,"flow_change/TTPair_VS_TTPco2_bg.jpeg", width=800, height=800)
dev.off()

# all_regression wt ntx     NOT CENTRED
wt_ntx_flow_change_tab = flowChange(data_matched,"wt&ntx")
# all_regression wt tx     NOT CENTRED
wt_tx_flow_change_tab = flowChange(data_matched,"wt&tx")
# all_regression tg ntx     NOT CENTRED
tg_ntx_flow_change_tab = flowChange(data_matched,"tg&ntx")
# all_regression tg tx     NOT CENTRED
tg_tx_flow_change_tab = flowChange(data_matched,"tg&tx")

p_wtVSwtLNAME <- sigP(data_matched,"wt&ntx","wt&tx")
p_wtVStg <- sigP(data_matched,"wt&ntx","tg&ntx")
p_tgVStgLNAME <- sigP(data_matched,"tg&ntx","tg&tx")

# art_regression wt ntx     NOT CENTRED
wt_ntx_art_flow_change_tab = flowChange(subset(data_matched,type=="a"),"wt&ntx")
# art_regression wt tx     NOT CENTRED
wt_tx_art_flow_change_tab = flowChange(subset(data_matched,type=="a"),"wt&tx")
# art_regression tg ntx     NOT CENTRED
tg_ntx_art_flow_change_tab = flowChange(subset(data_matched,type=="a"),"tg&ntx")
# art_regression tg tx     NOT CENTRED
tg_tx_art_flow_change_tab = flowChange(subset(data_matched,type=="a"),"tg&tx")

# cap_regression wt ntx     NOT CENTRED
wt_ntx_cap_flow_change_tab = flowChange(subset(data_matched,type=="c"),"wt&ntx")
# cap_regression wt tx     NOT CENTRED
wt_tx_cap_flow_change_tab = flowChange(subset(data_matched,type=="c"),"wt&tx")
# cap_regression tg ntx     NOT CENTRED
tg_ntx_cap_flow_change_tab = flowChange(subset(data_matched,type=="c"),"tg&ntx")
# cap_regression tg tx     NOT CENTRED
tg_tx_cap_flow_change_tab = flowChange(subset(data_matched,type=="c"),"tg&tx")

# ven_regression wt ntx     NOT CENTRED
wt_ntx_ven_flow_change_tab = flowChange(subset(data_matched,type=="v"),"wt&ntx")
# ven_regression wt tx     NOT CENTRED
wt_tx_ven_flow_change_tab = flowChange(subset(data_matched,type=="v"),"wt&tx")
# ven_regression tg ntx     NOT CENTRED
tg_ntx_ven_flow_change_tab = flowChange(subset(data_matched,type=="v"),"tg&ntx")
# ven_regression tg tx     NOT CENTRED
tg_tx_ven_flow_change_tab = flowChange(subset(data_matched,type=="v"),"tg&tx")

# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_flow_change = (tg_ntx_pca$coefficient[2]-wt_ntx_pca$coefficient[2])/
  (sqrt(tg_ntx_slp_se^2+wt_ntx_slp_se^2))
p_ntx_flow_change=2*pnorm(-abs(z_ntx_flow_change))
# Combined plot of all vessels of nTg and all vessels of Tg  
combined_plot <- ggplotRegressionPCAcombined(tg_ntx_ac,wt_ntx_ac,'darkorange3','olivedrab4',5,12,5.5,12)
combined_plot

bnd.df <- data.frame(matrix(ncol=4,nrow=11))
colnames(rfc.df) <- c("type","gn","rfc","se")
rfc.df[,3]<-c(wt_ntx_flow_change,tg_ntx_flow_change,wt_ntx_art_flow_change,tg_ntx_art_flow_change,
              wt_ntx_cap_flow_change,tg_ntx_cap_flow_change,wt_ntx_ven_flow_change,tg_ntx_ven_flow_change)
rfc.df[,4]<-c(wt_ntx_flow_se,tg_ntx_flow_se,wt_ntx_art_flow_se,tg_ntx_art_flow_se,
              wt_ntx_cap_flow_se,tg_ntx_cap_flow_se,wt_ntx_ven_flow_se,tg_ntx_ven_flow_se)
rfc.df[,4]<-(rfc.df[,4]/rfc.df[,3])*100     # se of % change in flow from co2 challenge
rfc.df[,3]<-(rfc.df[,3]-1)*100              # % change in flow from co2 challenge
rfc.df[,1] <- c("All Vessels","All Vessels","Arterioles","Arterioles",
                "Capillaries","Capillaries","Venules","Venules")
rfc.df[,2] <- c("nTg","Tg","nTg","Tg","nTg","Tg","nTg","Tg")
x.seq <- c(2,4,6.5,7.5,9.5,10.5,12.5,13.5)

comb <- combined_plot+
  labs(x=expression(paste(TTP[air], " (s)")), y=expression(paste(TTP[CO[2]], " (s)")))+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed",size=0.8)+
  #   geom_abline(slope=0,intercept=4,colour="blue",linetype="dashed",size=0.8)+
  geom_ribbon(aes(ymin=cl,ymax=cu),fill="darkgrey",alpha="0.5")+
  geom_ribbon(aes(ymin=cu,ymax=cl),fill="darkgrey",alpha="0.5")+
  #   scale_color_manual(name="Genotype",values=c("darkorange3","olivedrab4"),
  #                     labels=c("TgAD","nTg"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+theme_gray()+
  theme(axis.text=element_text(size=30),text=element_text(size=30),
        legend.position="top")+
  xlim(4,13)+ylim(4,13)
comb
dev.copy(jpeg,"flow_change/TTPair_VS_TTPco2.jpeg", width=800, height=800)
dev.off()

# Arterioles_regression wt ntx
airdata_wt_ntx_art_mtchd<-subset(airdata_wt_ntx_mtchd,type=="a")
cdata_wt_ntx_art_mtchd<-subset(cdata_wt_ntx_mtchd,type=="a")
wt_ntx_art_ac <-data.frame(airdata_wt_ntx_art_mtchd$pt,cdata_wt_ntx_art_mtchd$pt)
colnames(wt_ntx_art_ac) <- c("air","co2")
wt_ntx_art_pca <- deming(wt_ntx_art_ac[,2]~wt_ntx_art_ac[,1])
wt_ntx_art_flow_change <- 1/wt_ntx_art_pca$coefficients[2]
wt_ntx_art_ci_int <- (wt_ntx_art_pca$ci[2,2]-wt_ntx_art_pca$ci[2,1])/2
wt_ntx_art_flow_change_ci <- wt_ntx_art_ci_int*(wt_ntx_art_flow_change)^2
wt_ntx_art_slp_se <- sqrt(diag(wt_ntx_art_pca$variance))[2]                     # se for slope
wt_ntx_art_flow_se <- wt_ntx_art_slp_se/wt_ntx_art_pca$coefficients[2]^2        # se for flow change
wt_ntx_art_slp_int_se <- sqrt(diag(wt_ntx_art_pca$variance))[1]                 # se for original slope's intercept
wt_ntx_art_flow_int_se <- wt_ntx_art_slp_int_se/wt_ntx_art_pca$coefficients[1]^2# se for flow change
wt_ntx_art_ac$se_u <- (wt_ntx_art_pca$coefficients[2]+wt_ntx_art_slp_se)*wt_ntx_art_ac$air+wt_ntx_art_pca$coeff[1]-wt_ntx_slp_int_se
wt_ntx_art_ac$se_l <- (wt_ntx_art_pca$coefficients[2]-wt_ntx_art_slp_se)*wt_ntx_art_ac$air+wt_ntx_art_pca$coeff[1]+wt_ntx_slp_int_se
wt_ntx_art_ac$cu <- wt_ntx_art_pca$ci[2,2]*wt_ntx_art_ac$air+wt_ntx_art_pca$ci[1,1]
wt_ntx_art_ac$cl <- wt_ntx_art_pca$ci[2,1]*wt_ntx_art_ac$air+wt_ntx_art_pca$ci[1,2]
wt_ntx_art_plot <- ggplotRegressionPCA(wt_ntx_art_ac,5.5,12)
aa <- wt_ntx_art_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("nTg - Arterioles")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_nTG_nTX_Arterioles.jpeg")
dev.off()
wt_ntx_art_flow_change_tab = paste(round(wt_ntx_art_flow_change,digits=3)," +/- ",round(wt_ntx_art_flow_se,digits=3))
# Arterioles_regression tg ntx
airdata_tg_ntx_art_mtchd<-subset(airdata_tg_ntx_mtchd,type=="a")
cdata_tg_ntx_art_mtchd<-subset(cdata_tg_ntx_mtchd,type=="a")
tg_ntx_art_ac <-data.frame(airdata_tg_ntx_art_mtchd$pt,cdata_tg_ntx_art_mtchd$pt)
colnames(tg_ntx_art_ac) <- c("air","co2")
tg_ntx_art_pca <- deming(tg_ntx_art_ac[,2]~tg_ntx_art_ac[,1])
tg_ntx_art_flow_change <- 1/tg_ntx_art_pca$coefficients[2]
tg_ntx_art_ci_int <- (tg_ntx_art_pca$ci[2,2]-tg_ntx_art_pca$ci[2,1])/2
tg_ntx_art_flow_change_ci <- tg_ntx_art_ci_int*(tg_ntx_art_flow_change)^2
tg_ntx_art_slp_se <- sqrt(diag(tg_ntx_art_pca$variance))[2]                     # se for slope
tg_ntx_art_flow_se <- tg_ntx_art_slp_se/tg_ntx_art_pca$coefficients[2]^2        # se for flow change
tg_ntx_art_slp_int_se <- sqrt(diag(tg_ntx_art_pca$variance))[1]                 # se for original slope's intercept
tg_ntx_art_flow_int_se <- tg_ntx_art_slp_int_se/tg_ntx_art_pca$coefficients[1]^2# se for flow change
tg_ntx_art_ac$se_u <- (tg_ntx_art_pca$coefficients[2]+tg_ntx_art_slp_se)*tg_ntx_art_ac$air+tg_ntx_art_pca$coeff[1]-tg_ntx_slp_int_se
tg_ntx_art_ac$se_l <- (tg_ntx_art_pca$coefficients[2]-tg_ntx_art_slp_se)*tg_ntx_art_ac$air+tg_ntx_art_pca$coeff[1]+tg_ntx_slp_int_se
tg_ntx_art_ac$cu <- tg_ntx_art_pca$ci[2,2]*tg_ntx_art_ac$air+tg_ntx_art_pca$ci[1,1]
tg_ntx_art_ac$cl <- tg_ntx_art_pca$ci[2,1]*tg_ntx_art_ac$air+tg_ntx_art_pca$ci[1,2]
tg_ntx_art_plot <- ggplotRegressionPCA(tg_ntx_art_ac,5,12)
aa <- tg_ntx_art_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("Tg - Arterioles")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_TG_nTX_Arterioles.jpeg")
dev.off()
tg_ntx_art_flow_change_tab = paste(round(tg_ntx_art_flow_change,digits=3)," +/- ",round(tg_ntx_art_flow_se,digits=3))
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
wt_ntx_cap_flow_change_ci <- wt_ntx_cap_ci_int*(wt_ntx_cap_flow_change)^2
wt_ntx_cap_slp_se <- sqrt(diag(wt_ntx_cap_pca$variance))[2]                     # se for slope
wt_ntx_cap_flow_se <- wt_ntx_cap_slp_se/wt_ntx_cap_pca$coefficients[2]^2        # se for flow change
wt_ntx_cap_slp_int_se <- sqrt(diag(wt_ntx_cap_pca$variance))[1]                 # se for original slope's intercept
wt_ntx_cap_flow_int_se <- wt_ntx_cap_slp_int_se/wt_ntx_cap_pca$coefficients[1]^2# se for flow change
wt_ntx_cap_ac$se_u <- (wt_ntx_cap_pca$coefficients[2]+wt_ntx_cap_slp_se)*wt_ntx_cap_ac$air+wt_ntx_cap_pca$coeff[1]-wt_ntx_slp_int_se
wt_ntx_cap_ac$se_l <- (wt_ntx_cap_pca$coefficients[2]-wt_ntx_cap_slp_se)*wt_ntx_cap_ac$air+wt_ntx_cap_pca$coeff[1]+wt_ntx_slp_int_se
wt_ntx_cap_ac$cu <- wt_ntx_cap_pca$ci[2,2]*wt_ntx_cap_ac$air+wt_ntx_cap_pca$ci[1,1]
wt_ntx_cap_ac$cl <- wt_ntx_cap_pca$ci[2,1]*wt_ntx_cap_ac$air+wt_ntx_cap_pca$ci[1,2]
wt_ntx_cap_plot <- ggplotRegressionPCA(wt_ntx_cap_ac,5.5,12)
aa <- wt_ntx_cap_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("nTg - Capillaries")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_nTG_nTX_Capillaries.jpeg")
dev.off()
wt_ntx_cap_flow_change_tab = paste(round(wt_ntx_cap_flow_change,digits=3)," +/- ",round(wt_ntx_cap_flow_se,digits=3))
# Capillaries_regression tg ntx
airdata_tg_ntx_cap_mtchd<-subset(airdata_tg_ntx_mtchd,type=="c")
cdata_tg_ntx_cap_mtchd<-subset(cdata_tg_ntx_mtchd,type=="c")
tg_ntx_cap_ac <-data.frame(airdata_tg_ntx_cap_mtchd$pt,cdata_tg_ntx_cap_mtchd$pt)
colnames(tg_ntx_cap_ac) <- c("air","co2")
tg_ntx_cap_pca <- deming(tg_ntx_cap_ac[,2]~tg_ntx_cap_ac[,1])
tg_ntx_cap_flow_change <- 1/tg_ntx_cap_pca$coefficients[2]
tg_ntx_cap_ci_int <- (tg_ntx_cap_pca$ci[2,2]-tg_ntx_cap_pca$ci[2,1])/2
tg_ntx_cap_flow_change_ci <- tg_ntx_cap_ci_int*(tg_ntx_cap_flow_change)^2
tg_ntx_cap_slp_se <- sqrt(diag(tg_ntx_cap_pca$variance))[2]                     # se for slope
tg_ntx_cap_flow_se <- tg_ntx_cap_slp_se/tg_ntx_cap_pca$coefficients[2]^2        # se for flow change
tg_ntx_cap_slp_int_se <- sqrt(diag(tg_ntx_cap_pca$variance))[1]                 # se for original slope's intercept
tg_ntx_cap_flow_int_se <- tg_ntx_cap_slp_int_se/tg_ntx_cap_pca$coefficients[1]^2# se for flow change
tg_ntx_cap_ac$se_u <- (tg_ntx_cap_pca$coefficients[2]+tg_ntx_cap_slp_se)*tg_ntx_cap_ac$air+tg_ntx_cap_pca$coeff[1]-tg_ntx_slp_int_se
tg_ntx_cap_ac$se_l <- (tg_ntx_cap_pca$coefficients[2]-tg_ntx_cap_slp_se)*tg_ntx_cap_ac$air+tg_ntx_cap_pca$coeff[1]+tg_ntx_slp_int_se
tg_ntx_cap_ac$cu <- tg_ntx_cap_pca$ci[2,2]*tg_ntx_cap_ac$air+tg_ntx_cap_pca$ci[1,1]
tg_ntx_cap_ac$cl <- tg_ntx_cap_pca$ci[2,1]*tg_ntx_cap_ac$air+tg_ntx_cap_pca$ci[1,2]
tg_ntx_cap_plot <- ggplotRegressionPCA(tg_ntx_cap_ac,5,12)
aa <- tg_ntx_cap_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("Tg - Capillaries")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_TG_nTX_Capillaries.jpeg")
dev.off()
tg_ntx_cap_flow_change_tab = paste(round(tg_ntx_cap_flow_change,digits=3)," +/- ",round(tg_ntx_cap_flow_se,digits=3))
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
wt_ntx_ven_flow_change_ci <- wt_ntx_ven_ci_int*(wt_ntx_ven_flow_change)^2
wt_ntx_ven_slp_se <- sqrt(diag(wt_ntx_ven_pca$variance))[2]                     # se for slope
wt_ntx_ven_flow_se <- wt_ntx_ven_slp_se/wt_ntx_ven_pca$coefficients[2]^2        # se for flow change
wt_ntx_ven_slp_int_se <- sqrt(diag(wt_ntx_ven_pca$variance))[1]                 # se for original slope's intercept
wt_ntx_ven_flow_int_se <- wt_ntx_ven_slp_int_se/wt_ntx_ven_pca$coefficients[1]^2# se for flow change
wt_ntx_ven_ac$se_u <- (wt_ntx_ven_pca$coefficients[2]+wt_ntx_ven_slp_se)*wt_ntx_ven_ac$air+wt_ntx_ven_pca$coeff[1]-wt_ntx_slp_int_se
wt_ntx_ven_ac$se_l <- (wt_ntx_ven_pca$coefficients[2]-wt_ntx_ven_slp_se)*wt_ntx_ven_ac$air+wt_ntx_ven_pca$coeff[1]+wt_ntx_slp_int_se
wt_ntx_ven_ac$cu <- wt_ntx_ven_pca$ci[2,2]*wt_ntx_ven_ac$air+wt_ntx_ven_pca$ci[1,1]
wt_ntx_ven_ac$cl <- wt_ntx_ven_pca$ci[2,1]*wt_ntx_ven_ac$air+wt_ntx_ven_pca$ci[1,2]
wt_ntx_ven_plot <- ggplotRegressionPCA(wt_ntx_ven_ac,5.5,12)
aa <- wt_ntx_ven_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("nTg - Venules")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_nTG_nTX_Venules.jpeg")
dev.off()
wt_ntx_ven_flow_change_tab = paste(round(wt_ntx_ven_flow_change,digits=3)," +/- ",round(wt_ntx_ven_flow_se,digits=3))
# Venules_regression tg ntx
airdata_tg_ntx_ven_mtchd<-subset(airdata_tg_ntx_mtchd,type=="v")
cdata_tg_ntx_ven_mtchd<-subset(cdata_tg_ntx_mtchd,type=="v")
tg_ntx_ven_ac <-data.frame(airdata_tg_ntx_ven_mtchd$pt,cdata_tg_ntx_ven_mtchd$pt)
colnames(tg_ntx_ven_ac) <- c("air","co2")
tg_ntx_ven_pca <- deming(tg_ntx_ven_ac[,2]~tg_ntx_ven_ac[,1])
tg_ntx_ven_flow_change <- 1/tg_ntx_ven_pca$coefficients[2]
tg_ntx_ven_ci_int <- (tg_ntx_ven_pca$ci[2,2]-tg_ntx_ven_pca$ci[2,1])/2
tg_ntx_ven_flow_change_ci <- tg_ntx_ven_ci_int*(tg_ntx_ven_flow_change)^2
tg_ntx_ven_slp_se <- sqrt(diag(tg_ntx_ven_pca$variance))[2]                     # se for slope
tg_ntx_ven_flow_se <- tg_ntx_ven_slp_se/tg_ntx_ven_pca$coefficients[2]^2        # se for flow change
tg_ntx_ven_slp_int_se <- sqrt(diag(tg_ntx_ven_pca$variance))[1]                 # se for original slope's intercept
tg_ntx_ven_flow_int_se <- tg_ntx_ven_slp_int_se/tg_ntx_ven_pca$coefficients[1]^2# se for flow change
tg_ntx_ven_ac$se_u <- (tg_ntx_ven_pca$coefficients[2]+tg_ntx_ven_slp_se)*tg_ntx_ven_ac$air+tg_ntx_ven_pca$coeff[1]-tg_ntx_slp_int_se
tg_ntx_ven_ac$se_l <- (tg_ntx_ven_pca$coefficients[2]-tg_ntx_ven_slp_se)*tg_ntx_ven_ac$air+tg_ntx_ven_pca$coeff[1]+tg_ntx_slp_int_se
tg_ntx_ven_ac$cu <- tg_ntx_ven_pca$ci[2,2]*tg_ntx_ven_ac$air+tg_ntx_ven_pca$ci[1,1]
tg_ntx_ven_ac$cl <- tg_ntx_ven_pca$ci[2,1]*tg_ntx_ven_ac$air+tg_ntx_ven_pca$ci[1,2]
tg_ntx_ven_plot <- ggplotRegressionPCA(tg_ntx_ven_ac,5.5,12)
aa <- tg_ntx_ven_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("Tg - Venules")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_TG_nTX_Venules.jpeg")
dev.off()
tg_ntx_ven_flow_change_tab = paste(round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_se,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_ven_flow_change = (tg_ntx_ven_pca$coefficient[2]-wt_ntx_ven_pca$coefficient[2])/
  (sqrt(tg_ntx_ven_slp_se^2+wt_ntx_ven_slp_se^2))
p_ntx_ven_flow_change=2*pnorm(-abs(z_ntx_ven_flow_change))

wt_ntx_flow_tabs <- data.frame(c(wt_ntx_flow_change_tab,wt_ntx_art_flow_change_tab,wt_ntx_cap_flow_change_tab,wt_ntx_ven_flow_change_tab))
tg_ntx_flow_tabs <- data.frame(c(tg_ntx_flow_change_tab,tg_ntx_art_flow_change_tab,tg_ntx_cap_flow_change_tab,tg_ntx_ven_flow_change_tab))
p_values <- data.frame(c(p_ntx_flow_change,p_ntx_art_flow_change,p_ntx_cap_flow_change,p_ntx_ven_flow_change))
flow_change_tabss <- cbind(wt_ntx_flow_tabs,tg_ntx_flow_tabs,p_values)
colnames(flow_change_tabss) <- c("nTg Flow Change","Tg Flow Change","p-value")
rownames(flow_change_tabss) <- c("ALL","Arterioles","Capillaries","Venules")
flow_change_tabs

# % change in flow -> (RFC-1)*100 %
rfc.df <- data.frame(matrix(ncol=4,nrow=8))
colnames(rfc.df) <- c("type","gn","rfc","se")
rfc.df[,3]<-c(wt_ntx_flow_change,tg_ntx_flow_change,wt_ntx_art_flow_change,tg_ntx_art_flow_change,
              wt_ntx_cap_flow_change,tg_ntx_cap_flow_change,wt_ntx_ven_flow_change,tg_ntx_ven_flow_change)
rfc.df[,4]<-c(wt_ntx_flow_se,tg_ntx_flow_se,wt_ntx_art_flow_se,tg_ntx_art_flow_se,
              wt_ntx_cap_flow_se,tg_ntx_cap_flow_se,wt_ntx_ven_flow_se,tg_ntx_ven_flow_se)
rfc.df[,4]<-(rfc.df[,4]/rfc.df[,3])*100     # se of % change in flow from co2 challenge
rfc.df[,3]<-(rfc.df[,3]-1)*100              # % change in flow from co2 challenge
rfc.df[,1] <- c("All Vessels","All Vessels","Arterioles","Arterioles",
                "Capillaries","Capillaries","Venules","Venules")
rfc.df[,2] <- c("nTg","Tg","nTg","Tg","nTg","Tg","nTg","Tg")
x.seq <- c(2,4,6.5,7.5,9.5,10.5,12.5,13.5)
rfc.df$x <- x.seq
rfc.df$width[rfc.df$type=="All Vessels"] <- 2
rfc.df$width[rfc.df$type!="All Vessels"] <- 1
perc_change_temp<-ggplot(rfc.df, aes(x=x, y=rfc)) +
  geom_errorbar(data=subset(rfc.df,rfc>0),aes(fill=gn,ymin=rfc, ymax=rfc+se),
                size=.7,
                width=.4,
                position=position_dodge(0.9))+
  #   geom_errorbar(data=subset(rfc.df,rfc<0),aes(fill=gn,ymin=rfc, ymax=rfc-se),
  #                 size=.7,
  #                 width=.4,
  #                 position=position_dodge(0.9))+
  geom_bar(aes(fill=gn, width=width),position=position_dodge(), stat="identity",
           color="black",
           size=1)+ylim(-25,130)+xlab("")+ylab(expression(paste(Delta,"Flow [%]")))+
  scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","Tg"))+
  #   scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","Tg"))+
  scale_x_continuous(breaks=c(sum(x.seq[1:2])/2, sum(x.seq[3:4])/2, sum(x.seq[5:6])/2, sum(x.seq[7:8])/2),
                     label=c("All Vessels","Arterioles","Capillaries","Venules"))+theme_gray()+
  theme(axis.text=element_text(size=25),
        text=element_text(size=30),
        legend.text=element_text(size=30),
        legend.key.height=unit(2.5,"line"),
        legend.position="bottom")
perc_change_temp
# label1.df <- data.frame(x = c(sum(x.seq[1:2])/2,sum(x.seq[3:4])/2,sum(x.seq[5:6])/2),
#                         rfc = c(2.25,2.45,2.75))
# data1 <- data.frame(x = c(2,2,4,4), y = c(2,2.25,2.25,2),
#                     type=NA, delta_pt=NA)
# data2 <- data.frame(x = c(6.5,6.5,7.5,7.5), y = c(2.2,2.45,2.45,2.2),
#                     type=NA, delta_pt=NA)
# data3 <- data.frame(x = c(9.5,9.5,10.5,10.5), y = c(2.5,2.75,2.75,2.5),
#                     type=NA, delta_pt=NA)
label1.df <- data.frame(x = c(sum(x.seq[1:2])/2,sum(x.seq[3:4])/2,sum(x.seq[5:6])/2),
                        rfc = c(100,115,122))
data1 <- data.frame(x = c(2,2,4,4), y = c(90,100,100,90),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(6.5,6.5,7.5,7.5), y = c(105,115,115,105),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(9.5,9.5,10.5,10.5), y = c(112,122,122,112),
                    type=NA, delta_pt=NA)
perc_change<-perc_change_temp + geom_path(data = data1, aes(x = x, y = y),size=1.5)+
  geom_path(data = data2, aes(x = x, y = y),size=1.5)+
  geom_path(data = data3, aes(x = x, y = y),size=1.5)+
  geom_text(data=label1.df,label="*",size=18) + theme(legend.position="bottom")
perc_change
dev.copy(jpeg,"flow_change/perc_change_flow.jpeg", width=800, height=800)
dev.off()

# combining flow change plot ('comb') and percent change plot('perc_change')
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(perc_change)
perc_change <- perc_change + theme(legend.position="none")
blankPlot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()
top_row <- plot_grid(comb,blankPlot,perc_change, labels=c("A","","B"),label_size=25,ncol=3, rel_widths = c(1,0.1,1))
plot_grid(top_row, legend, ncol = 1, rel_heights = c(1, 0.1))
dev.copy(jpeg,"flow_change/airVSco2&perc_change1.jpeg", width=1600, height=800)
dev.off()
