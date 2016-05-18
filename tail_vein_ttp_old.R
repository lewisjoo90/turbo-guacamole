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
data<-as.data.frame(read.csv("bolus_AD_Rat_all_OLrej_overlap.csv",header=T))
#data_exc <- data[-c(15,16,39,40,41,42,45,46,51,52,57,58), ]

#consider these categorical variables
data$id <- factor(data$id)
# data$trt <- factor(data$trt)
data$vessel_number <- factor(data$vessel_number)
data$bolus_number <- factor(data$bolus_number)
data$state <- factor(data$state,levels=c("air","co2"))
data$depth <- factor(data$depth)
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
airdata_wt_ntx<-subset(data,data$state=="air"&data$gntrt=="wt&ntx")
airdata_tg_ntx<-subset(data,data$state=="air"&data$gntrt=="tg&ntx")
cdata_wt_ntx<-subset(data,data$state=="co2"&data$gntrt=="wt&ntx")
cdata_tg_ntx<-subset(data,data$state=="co2"&data$gntrt=="tg&ntx")
airdata_wt_tx<-subset(data,data$state=="air"&data$gntrt=="wt&tx")
airdata_tg_tx<-subset(data,data$state=="air"&data$gntrt=="tg&tx")
cdata_wt_tx<-subset(data,data$state=="co2"&data$gntrt=="wt&tx")
cdata_tg_tx<-subset(data,data$state=="co2"&data$gntrt=="tg&tx")

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

## break down the data by vessel type
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
sd_ntx_wt_air_a <- sd(airdata_ntx_wt[airdata_ntx_wt$type=="a",]$pt)
sd_ntx_wt_co2_a <- sd(cdata_ntx_wt[cdata_ntx_wt$type=="a",]$pt)
sd_ntx_wt_air_c <- sd(airdata_ntx_wt[airdata_ntx_wt$type=="c",]$pt)
sd_ntx_wt_co2_c <- sd(cdata_ntx_wt[cdata_ntx_wt$type=="c",]$pt)
sd_ntx_wt_air_v <- sd(airdata_ntx_wt[airdata_ntx_wt$type=="v",]$pt)
sd_ntx_wt_co2_v <- sd(cdata_ntx_wt[cdata_ntx_wt$type=="v",]$pt)
sd_ntx_tg_air_a <- sd(airdata_ntx_tg[airdata_ntx_tg$type=="a",]$pt)
sd_ntx_tg_co2_a <- sd(cdata_ntx_tg[cdata_ntx_tg$type=="a",]$pt)
sd_ntx_tg_air_c <- sd(airdata_ntx_tg[airdata_ntx_tg$type=="c",]$pt)
sd_ntx_tg_co2_c <- sd(cdata_ntx_tg[cdata_ntx_tg$type=="c",]$pt)
sd_ntx_tg_air_v <- sd(airdata_ntx_tg[airdata_ntx_tg$type=="v",]$pt)
sd_ntx_tg_co2_v <- sd(cdata_ntx_tg[cdata_ntx_tg$type=="v",]$pt)
mean_ntx_wt_air_a <- mean(airdata_ntx_wt[airdata_ntx_wt$type=="a",]$pt)
mean_ntx_wt_co2_a <- mean(cdata_ntx_wt[cdata_ntx_wt$type=="a",]$pt)
mean_ntx_wt_air_c <- mean(airdata_ntx_wt[airdata_ntx_wt$type=="c",]$pt)
mean_ntx_wt_co2_c <- mean(cdata_ntx_wt[cdata_ntx_wt$type=="c",]$pt)
mean_ntx_wt_air_v <- mean(airdata_ntx_wt[airdata_ntx_wt$type=="v",]$pt)
mean_ntx_wt_co2_v <- mean(cdata_ntx_wt[cdata_ntx_wt$type=="v",]$pt)
mean_ntx_tg_air_a <- mean(airdata_ntx_tg[airdata_ntx_tg$type=="a",]$pt)
mean_ntx_tg_co2_a <- mean(cdata_ntx_tg[cdata_ntx_tg$type=="a",]$pt)
mean_ntx_tg_air_c <- mean(airdata_ntx_tg[airdata_ntx_tg$type=="c",]$pt)
mean_ntx_tg_co2_c <- mean(cdata_ntx_tg[cdata_ntx_tg$type=="c",]$pt)
mean_ntx_tg_air_v <- mean(airdata_ntx_tg[airdata_ntx_tg$type=="v",]$pt)
mean_ntx_tg_co2_v <- mean(cdata_ntx_tg[cdata_ntx_tg$type=="v",]$pt)
cov_ntx_wt_air_a <- sd_ntx_wt_air_a/mean_ntx_wt_air_a
cov_ntx_wt_co2_a <- sd_ntx_wt_co2_a/mean_ntx_wt_co2_a
cov_ntx_wt_air_c <- sd_ntx_wt_air_c/mean_ntx_wt_air_c
cov_ntx_wt_co2_c <- sd_ntx_wt_co2_c/mean_ntx_wt_co2_c
cov_ntx_wt_air_v <- sd_ntx_wt_air_v/mean_ntx_wt_air_v
cov_ntx_wt_co2_v <- sd_ntx_wt_co2_v/mean_ntx_wt_co2_v
cov_ntx_tg_air_a <- sd_ntx_tg_air_a/mean_ntx_tg_air_a
cov_ntx_tg_co2_a <- sd_ntx_tg_co2_a/mean_ntx_tg_co2_a
cov_ntx_tg_air_c <- sd_ntx_tg_air_c/mean_ntx_tg_air_c
cov_ntx_tg_co2_c <- sd_ntx_tg_co2_c/mean_ntx_tg_co2_c
cov_ntx_tg_air_v <- sd_ntx_tg_air_v/mean_ntx_tg_air_v
cov_ntx_tg_co2_v <- sd_ntx_tg_co2_v/mean_ntx_tg_co2_v
sd_wt <- c(sd_ntx_wt_air_a,sd_ntx_wt_co2_a,sd_ntx_wt_air_c,sd_ntx_wt_co2_c,
           sd_ntx_wt_air_v,sd_ntx_wt_co2_v)
sd_tg <- c(sd_ntx_tg_air_a,sd_ntx_tg_co2_a,sd_ntx_tg_air_c,sd_ntx_tg_co2_c,
           sd_ntx_tg_air_v,sd_ntx_tg_co2_v)
sd_all <- t(data.frame(sd_wt,sd_tg))
cov_wt <- c(cov_ntx_wt_air_a,cov_ntx_wt_co2_a,cov_ntx_wt_air_c,cov_ntx_wt_co2_c,
           cov_ntx_wt_air_v,cov_ntx_wt_co2_v)
cov_tg <- c(cov_ntx_tg_air_a,cov_ntx_tg_co2_a,cov_ntx_tg_air_c,cov_ntx_tg_co2_c,
           cov_ntx_tg_air_v,cov_ntx_tg_co2_v)
cov_all <- t(data.frame(cov_wt,cov_tg))
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="a",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="c",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="v",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="v",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="a",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="c",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="v",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="v",]$pt,alternative=c("greater"))

FC_60sec <- c(paste(round(pct50_60sec_flow_change,digits=3)," +/- ",round(pct50_60sec_flow_change_std,digits=3)),
              paste(round(pct75_60sec_flow_change,digits=3)," +/- ",round(pct75_60sec_flow_change_std,digits=3)),
              paste(round(pct100_60sec_flow_change,digits=3)," +/- ",round(pct100_60sec_flow_change_std,digits=3)))
FC_60sec_cap <- c(paste(round(pct50_60sec_cap_flow_change,digits=3)," +/- ",round(pct50_60sec_cap_flow_change_std,digits=3)),
                  paste(round(pct75_60sec_cap_flow_change,digits=3)," +/- ",round(pct75_60sec_cap_flow_change_std,digits=3)),
                  paste(round(pct100_60sec_cap_flow_change,digits=3)," +/- ",round(pct100_60sec_cap_flow_change_std,digits=3)))
FC_60sec_pv <- c(paste(round(pct50_60sec_pv_flow_change,digits=3)," +/- ",round(pct50_60sec_pv_flow_change_std,digits=3)),
                 paste(round(pct75_60sec_pv_flow_change,digits=3)," +/- ",round(pct75_60sec_pv_flow_change_std,digits=3)),
                 paste(round(pct100_60sec_pv_flow_change,digits=3)," +/- ",round(pct100_60sec_pv_flow_change_std,digits=3)))
FC_60sec_df <-data.frame(FC_60sec,FC_60sec_cap,FC_60sec_pv)
setattr(FC_60sec_df,"row.names",c("5% CO2 - 60sec","7.5% CO2 - 60sec","10% CO2 - 60sec"))

FC_60sec_df.table <- xtable(FC_60sec_df[1:3,])
print.xtable(FC_60sec_df.table,floating=FALSE,type="latex",file="aaa.tex")





# barplot change in ttp by genotype
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
# TTP box
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
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
aa = ggplot(wt_ntx_ac,aes(pt,colour=state))+geom_freqpoly(binwidth=0.5)
aa

# all_regression wt ntx
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
df <- wt_ntx_ac
colnames(df) <- c("air","co2")
wt_ntx_pca <- deming(df[,2]~df[,1])
line_cu <- list(slope=wt_ntx_pca$ci[2,2],intercept=wt_ntx_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_pca$ci[2,1],intercept=wt_ntx_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
df$cu <- wt_ntx_pca$ci[2,2]*df$air+wt_ntx_pca$ci[1,1]
df$cl <- wt_ntx_pca$ci[2,1]*df$air+wt_ntx_pca$ci[1,2]
cu_coord <- wt_ntx_pca$ci[2,2]*x_coord+wt_ntx_pca$ci[1,1]
cl_coord <- wt_ntx_pca$ci[2,1]*x_coord+wt_ntx_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_pca$ci[1,2]-wt_ntx_pca$ci[1,1])/(wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1]),
               y=(wt_ntx_pca$ci[2,1]*wt_ntx_pca$ci[1,1]-wt_ntx_pca$ci[2,2]*wt_ntx_pca$ci[1,2])/(wt_ntx_pca$ci[2,1]-wt_ntx_pca$ci[2,2])))
ggplot(df,aes_string(names(df)[1],names(df)[2]))+geom_point()+
  stat_smooth(data=subset(df,air>5.5&air<12),
              method = f, se = TRUE, colour = 'red',size=1.5, formula=y~x)+
  geom_line(aes(x=air,y=cu),colour="darkgrey",alpha="0.5")+
  geom_line(aes(x=air,y=cl),colour="darkgrey",alpha="0.5")+
  geom_ribbon(aes(ymin=cl,ymax=cu),fill="darkgrey",alpha="0.5")+
  geom_ribbon(aes(ymin=cu,ymax=cl),fill="darkgrey",alpha="0.5")+
  xlim(4,12.5)+ylim(4,12.5)
  geom_abline()

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
ggplotRegressionPCA <- function (df) {
  
  require(ggplot2)
  pca <- prcomp(df)$rotation
  pca_slp <- pca[2,1]/pca[1,1]
  
  ggplot(df, aes_string(x = names(df)[1], y = names(df)[2])) + 
    geom_point() +
    geom_abline(slope=pca_slp,colour="red") +
    labs(title = paste(" Slope =",signif(pca_slp, 2)))
}
ggplotRegressionPCA <- function (df) {
  
  require(ggplot2)
  pca <- prcomp(df)$rotation
  pca_slp <- pca[2,1]/pca[1,1]
  
  ggplot(df, aes_string(x = names(df)[1], y = names(df)[2])) + 
    geom_point() +
    stat_smooth(method="deming",colour="red") +
    labs(title = paste(" Slope =",signif(pca_slp, 2)))
}

# how to manually add ci in the ggplot: http://www.r-bloggers.com/shading-between-two-lines-ggplot/

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
# all_regression tg ntx
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

