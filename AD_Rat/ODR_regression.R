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
airdata_wt_ntx<-subset(data,data$state=="air"&data$gntrt=="wt&ntx")
airdata_tg_ntx<-subset(data,data$state=="air"&data$gntrt=="tg&ntx")
cdata_wt_ntx<-subset(data,data$state=="co2"&data$gntrt=="wt&ntx")
cdata_tg_ntx<-subset(data,data$state=="co2"&data$gntrt=="tg&ntx")
airdata_wt_tx<-subset(data,data$state=="air"&data$gntrt=="wt&tx")
airdata_tg_tx<-subset(data,data$state=="air"&data$gntrt=="tg&tx")
cdata_wt_tx<-subset(data,data$state=="co2"&data$gntrt=="wt&tx")
cdata_tg_tx<-subset(data,data$state=="co2"&data$gntrt=="tg&tx")

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
aa <- wt_ntx_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y="TTP during CO2 breathing [s]")+labs(x="TTP during air breathing [s]")+ggtitle("nTg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))
aa
dev.copy(jpeg,"FC_pca_nTG_nTX_ALL.jpeg")
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
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
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
dev.copy(jpeg,"FC_pca_TG_nTX_ALL.jpeg")
dev.off()
tg_ntx_flow_change_tab = paste(round(tg_ntx_flow_change,digits=3)," +/- ",round(tg_ntx_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
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
dev.copy(jpeg,"FC_pca_nTG_nTX_Arterioles.jpeg")
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
dev.copy(jpeg,"FC_pca_TG_nTX_Arterioles.jpeg")
dev.off()
tg_ntx_art_flow_change_tab = paste(round(tg_ntx_art_flow_change,digits=3)," +/- ",round(tg_ntx_art_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
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
dev.copy(jpeg,"FC_pca_nTG_nTX_Capillaries.jpeg")
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
dev.copy(jpeg,"FC_pca_TG_nTX_Capillaries.jpeg")
dev.off()
tg_ntx_cap_flow_change_tab = paste(round(tg_ntx_cap_flow_change,digits=3)," +/- ",round(tg_ntx_cap_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
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
dev.copy(jpeg,"FC_pca_nTG_nTX_Venules.jpeg")
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
dev.copy(jpeg,"FC_pca_TG_nTX_Venules.jpeg")
dev.off()
tg_ntx_ven_flow_change_tab = paste(round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_change_ci,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
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