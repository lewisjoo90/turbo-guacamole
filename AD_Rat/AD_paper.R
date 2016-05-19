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
library(cowplot)
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
ggplotRegressionPCAcombined <- function (df1,df2,colorVar1,colorVar2,x_left1,x_right1,x_left2,x_right2) {
  require(ggplot2)
  df3 <- data.frame(gn=c(rep("Tg",nrow(df1)), rep("nTg",nrow(df2))))
  df3 <- cbind(df3,rbind(df1,df2))
  ggplot(df3,aes_string(names(df3)[2],names(df3)[3]))+geom_point(aes(color=gn),size=2,alpha=0.4)+
    geom_line(data=df3,aes(x=air,y=cu,group=gn),colour="darkgrey",alpha="0.5")+
    geom_line(data=df3,aes(x=air,y=cl,group=gn),colour="darkgrey",alpha="0.5")+
    geom_ribbon(data=df3,aes(ymin=cl,ymax=cu,group=gn),fill="darkgrey",alpha=0.6)+
    geom_ribbon(data=df3,aes(ymin=cu,ymax=cl,group=gn),fill="darkgrey",alpha=0.6)+
    scale_color_manual(name="Genotype",values=c(colorVar2,colorVar1),
                       labels=c("nTg","Tg"))+
#     geom_point(data=df2,color=colorVar2,size=2,alpha=0.4)+
#     geom_line(data=df2,aes(x=air,y=cu),colour="darkgrey",alpha="0.5")+
#     geom_line(data=df2,aes(x=air,y=cl),colour="darkgrey",alpha="0.5")+
#     geom_ribbon(data=df2,aes(ymin=cl,ymax=cu),fill="darkgrey",alpha=0.4)+
#     geom_ribbon(data=df2,aes(ymin=cu,ymax=cl),fill="darkgrey",alpha=0.4)+
    stat_smooth(data=subset(df1,air>x_left1&air<x_right1),
                method = f, se = TRUE, colour = colorVar1, size=1.5, formula=y~x)+
    stat_smooth(data=subset(df2,air>x_left2&air<x_right2),
                method = f, se = TRUE, colour = colorVar2, size=1.5, formula=y~x)
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

# immunoblots
imm<-as.data.frame(read.csv("immuno.csv",header=T))
summary(lm(value ~ gn ,data = subset(imm,stain=="pdgfrb")))

# pdg <- ggplot(subset(imm,stain=="pdgfrb"),aes(x=gn,y=value,color=stain))+ylim(0,300)+labs(x="",y="")+
#   geom_jitter(position=position_jitter(w=0.1,h=0.1))
#   layer(data=data_wt_ntx_summ,mapping=aes(x=state,y=pt),geom="point",size=5,color="red")+
#   facet_grid(. ~ type,scales="free",space="free",labeller=vessel_labeller)+
#   scale_x_discrete(labels=c("Air","CO2"))+
#   ggtitle(paste("nTg"))+
#   theme(axis.text=element_text(size=16,colour="black"),
#         text=element_text(size=24),
#         legend.text=element_text(size=19))
# pdg

data_summ <- summarySE(imm, measurevar="value", groupvars=c("stain","gn"))
data_summ
x.seq <- c(1,2,1,2,1,2,1,2,1,2,1,2)
data_summ$x <- x.seq
data_summ

desmin<-ggplot(data=subset(data_summ,stain=="Desmin"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,500)+xlab("")+
  ylab(expression(atop("Desmin immunoblot density", paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle("Desmin")+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+
  theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
desmin
label1.df <- data.frame(x = c(x.seq[2]),
                        value = c(450))
Desmin<-desmin+
#   geom_path(data = data1, aes(x = x, y = y),size=1.5)+
#   geom_path(data = data2, aes(x = x, y = y),size=1.5)+
#   geom_path(data = data3, aes(x = x, y = y),size=1.5)+
#   geom_path(data = data4, aes(x = x, y = y),size=1.5)+
  geom_text(data=label1.df,label="*",size=18)
#   geom_text(data=label2.df,label="*",size=18)
Desmin
asma<-ggplot(data=subset(data_summ,stain=="a-sma"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,200)+xlab("")+
  ylab(expression(atop(paste(alpha,"-SMA immunoblot density"), paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle(expression(paste(alpha,"-SMA")))+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
asma
label1.df <- data.frame(x = c(x.seq[2]),
                        value = c(178))
Asma<-asma +
  #   geom_path(data = data1, aes(x = x, y = y),size=1.5)+
  #   geom_path(data = data2, aes(x = x, y = y),size=1.5)+
  #   geom_path(data = data3, aes(x = x, y = y),size=1.5)+
  #   geom_path(data = data4, aes(x = x, y = y),size=1.5)+
  geom_text(data=label1.df,label="*",size=18)
#   geom_text(data=label2.df,label="*",size=18)
Asma
pdg<-ggplot(data=subset(data_summ,stain=="pdgfrb"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,250)+xlab("")+
  ylab(expression(atop(paste("PDGFR",beta," immunoblot density"), paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle(expression(paste("PDGFR",beta)))+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
pdg
Pdg<-pdg
Pdg
ng2<-ggplot(data=subset(data_summ,stain=="ng2"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,150)+xlab("")+
  ylab(expression(atop(paste("NG2 immunoblot density"), paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle(expression(paste("NG2")))+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
ng2
Ng2<-ng2
Ng2
occ<-ggplot(data=subset(data_summ,stain=="occludin"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,150)+xlab("")+
  ylab(expression(atop(paste("Occludin immunoblot density"), paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle(expression(paste("Occludin")))+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
occ
Occ<-occ
Occ
gfa<-ggplot(data=subset(data_summ,stain=="gfap"),aes(x=x, y=value)) +
  geom_bar(aes(fill=gn),position=position_dodge(), stat="identity",
           color="black",
           size=1.2,width=0.6)+ylim(0,150)+xlab("")+
  ylab(expression(atop(paste("GFAP immunoblot density"), paste("normalized to GAPDH (% nTg)"))))+
  geom_errorbar(aes(fill=gn,ymin=value, ymax=value+se),
                size=.9,
                width=.4,
                position=position_dodge(0.9))+ggtitle(expression(paste("GFAP")))+
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","TgAD"))+
  scale_x_continuous(breaks=c(x.seq[1], x.seq[2]),label=c("nTg","Tg"))+theme_gray()+
  theme(axis.text=element_text(size=28),
        text=element_text(size=25),
        legend.text=element_text(size=20),
        legend.key.height=unit(2.5,"line"),
        legend.position="none")
gfa
Gfa<-gfa
Gfa

library(png)
library(grid)
img <- readPNG("immunoblots.PNG")
g <- rasterGrob(img, interpolate=TRUE)
immuno <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
immuno
immuno_width = 0.4
immuno_height = 0.5
stain_width = 0.25
stain_height = 0.4
ggdraw()+
  draw_plot(immuno,0,0.5,immuno_width,immuno_height)+
  draw_plot(Desmin,0.43,0.55,stain_width,stain_height)+
  draw_plot(Asma,0.71,0.55,stain_width,stain_height)+
  draw_plot(Pdg,0,0.05,stain_width,stain_height)+
  draw_plot(Ng2,0.25,0.05,stain_width,stain_height)+
  draw_plot(Occ,0.5,0.05,stain_width,stain_height)+
  draw_plot(Gfa,0.75,0.05,stain_width,stain_height)+
  draw_plot_label(c("A","B","C","D","E","F","G"),c(0,0.43,0.71,0,0.25,0.5,0.75),c(1,1,1,0.5,0.5,0.5,0.5),size=40)
dev.copy(jpeg,"histo+immuno/FIGURE2-barstest.jpeg", width=1600, height=1600)
dev.off()

  # barplot change in ttp by genotype for ntx only - ALL bar twice the width of single vessel type 
airdata_ntx<-subset(data_matched,data_matched$state=="air"&data_matched$trt=="ntx")
cdata_ntx<-subset(data_matched,data_matched$state=="co2"&data_matched$trt=="ntx")
delta_ntx_pt <- airdata_ntx$pt - cdata_ntx$pt
delta_p_ntx_pt <- (airdata_ntx$pt - cdata_ntx$pt)/airdata_ntx$pt*100
airdata_ntx$delta_pt <- delta_ntx_pt
airdata_ntx$delta_p_pt <- delta_p_ntx_pt
airdata_ntx <- airdata_ntx[order(rev(airdata_ntx$gn),airdata_ntx$id,airdata_ntx$vessel_number),] #
airdata_ntx <- subset(airdata_ntx,id!=74)
cdata_ntx <- cdata_ntx[order(rev(cdata_ntx$gn),cdata_ntx$id,cdata_ntx$vessel_number),] #
cdata_ntx <- subset(cdata_ntx,id!=74)
airdata_ntx$gn <- factor(airdata_ntx$gn,c("wt","tg"))
airdata_ntx_all <- airdata_ntx
cdata_ntx_all <- cdata_ntx
# For delta_pt vs vessel type including All Vessels
tmp <- as.character(airdata_ntx_all$type)
tmp <- "all"
airdata_ntx_all$type <- factor(tmp)
airdata_ntx_all <- rbind(airdata_ntx,airdata_ntx_all)
airdata_ntx_all$type = with(airdata_ntx_all, factor(type, levels = rev(levels(type))))
# delta_TTP bar 
data_summary <- summarySE(airdata_ntx_all, measurevar="delta_pt", groupvars=c("type","gn"))
data_summ <- data_summary
data_summ$width[data_summ$type=="all"] <- 2
data_summ$width[data_summ$type!="all"] <- 1
temp <- data_summ[3:4,]
data_summ[3:4,] <- data_summ[7:8,]
data_summ[7:8,] <- temp
x.seq <- c(2,4,6.5,7.5,9.5,10.5,12.5,13.5)
data_summ$x <- x.seq
data_summ
aa<-ggplot(data_summ, aes(x=x, y=delta_pt)) +
  geom_bar(aes(fill=gn, width=width),position=position_dodge(), stat="identity",
           color="black",
           size=1)+ylim(0,2.5)+xlab("")+ylab(expression(paste(Delta,"TTP [s]")))+
  geom_errorbar(aes(fill=gn,ymin=delta_pt, ymax=delta_pt+se),
                size=.7,
                width=.4,
                position=position_dodge(0.9)) +
#   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","Tg"))+
  scale_x_continuous(breaks=c(sum(x.seq[1:2])/2, sum(x.seq[3:4])/2, sum(x.seq[5:6])/2, sum(x.seq[7:8])/2),
                     label=c("All Vessels","Arterioles","Capillaries","Venules"))+theme_gray()+
  labs(y=expression(atop("Vascular Reactivity", paste("[",TTP[air]," - ",TTP[CO[2]],", s]"))))+
  theme(axis.text=element_text(size=25),
        text=element_text(size=30),
        legend.text=element_text(size=30),
        legend.key.height=unit(2.5,"line"),
        legend.position="bottom")
aa
label1.df <- data.frame(x = c(sum(x.seq[1:2])/2,sum(x.seq[3:4])/2,sum(x.seq[7:8])/2),
                        delta_pt = c(1.5,1.6,1.9))
label2.df <- data.frame(x = c(sum(x.seq[5:6])/2),
                        delta_pt = c(1.4))
data1 <- data.frame(x = c(2,2,4,4), y = c(1.4, 1.5, 1.5, 1.4),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(6.5,6.5,7.5,7.5), y = c(1.5, 1.6, 1.6, 1.5),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(9.5,9.5,10.5,10.5), y = c(1.3,1.4,1.4,1.3),
                    type=NA, delta_pt=NA)
data4 <- data.frame(x = c(12.5,12.5,13.5,13.5), y = c(1.8,1.9,1.9,1.8),
                    type=NA, delta_pt=NA)
bb<-aa + geom_path(data = data1, aes(x = x, y = y),size=1.5)+
  geom_path(data = data2, aes(x = x, y = y),size=1.5)+
#   geom_path(data = data3, aes(x = x, y = y),size=1.5)+
  geom_path(data = data4, aes(x = x, y = y),size=1.5)+
  geom_text(data=label1.df,label="*",size=18)
#   geom_text(data=label2.df,label="*",size=18)
bb
dev.copy(jpeg,"delta_ttp/delta_ttp_all_wtVStg1.jpeg", width=800, height=800)
dev.off()
airdata_ntx_all_all <- subset(airdata_ntx_all,type=="all")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_all)
airdata_ntx_all_art <- subset(airdata_ntx_all,type=="a")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_art)
airdata_ntx_all_cap <- subset(airdata_ntx_all,type=="c")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_cap)
airdata_ntx_all_ven <- subset(airdata_ntx_all,type=="v")
wilcox.test(delta_pt~gn,data=airdata_ntx_all_ven)

# delta%_TTP bar 
data_summary <- summarySE(airdata_ntx_all, measurevar="delta_p_pt", groupvars=c("type","gn"))
data_summ <- data_summary
data_summ$width[data_summ$type=="all"] <- 2
data_summ$width[data_summ$type!="all"] <- 1
temp <- data_summ[3:4,]
data_summ[3:4,] <- data_summ[7:8,]
data_summ[7:8,] <- temp
x.seq <- c(2,4,6.5,7.5,9.5,10.5,12.5,13.5)
data_summ$x <- x.seq
data_summ
aa<-ggplot(data_summ, aes(x=x, y=delta_p_pt)) +
  geom_bar(aes(fill=gn, width=width),position=position_dodge(), stat="identity",
           color="black",
           size=1)+ylim(0,30)+xlab("")+
  geom_errorbar(aes(fill=gn,ymin=delta_p_pt, ymax=delta_p_pt+se),
                size=.7,
                width=.4,
                position=position_dodge(0.9)) +
  #   scale_fill_manual(name="Genotype",values=c("olivedrab4","darkorange3"),labels=c("nTg","TgAD"))+
  scale_fill_manual(name="Genotype",values=c("white","grey45"),labels=c("nTg","Tg"))+
  scale_x_continuous(breaks=c(sum(x.seq[1:2])/2, sum(x.seq[3:4])/2, sum(x.seq[5:6])/2, sum(x.seq[7:8])/2),
                     label=c("All Vessels","Arterioles","Capillaries","Venules"))+theme_gray()+
#   labs(y=expression(atop("Vascular Reactivity [%]", paste("[",TTP[air]," - ",TTP[CO[2]],", %]"))))+
  labs(y="Vascular Reactivity [%]")+
  theme(axis.text=element_text(size=30),
        axis.title.y=element_text(vjust=5),
        text=element_text(size=30),
        legend.position="bottom")
aa
label1.df <- data.frame(x = c(sum(x.seq[1:2])/2,sum(x.seq[3:4])/2,sum(x.seq[7:8])/2),
                        delta_p_pt = c(20,21,24))
data1 <- data.frame(x = c(2,2,4,4), y = c(19, 20, 20, 19),
                    type=NA, delta_p_pt=NA)
data2 <- data.frame(x = c(6.5,6.5,7.5,7.5), y = c(20, 21, 21, 20),
                    type=NA, delta_p_pt=NA)
data4 <- data.frame(x = c(12.5,12.5,13.5,13.5), y = c(23,24,24,23),
                    type=NA, delta_p_pt=NA)
perc_ttp<-aa + geom_path(data = data1, aes(x = x, y = y),size=1.5)+
  geom_path(data = data2, aes(x = x, y = y),size=1.5)+
  #   geom_path(data = data3, aes(x = x, y = y),size=1.5)+
  geom_path(data = data4, aes(x = x, y = y),size=1.5)+
  geom_text(data=label1.df,label="*",size=18)
#   geom_text(data=label2.df,label="*",size=18)
perc_ttp
dev.copy(jpeg,"delta_ttp/dp_ttp_all_wtVStg1.jpeg", width=800, height=800)
dev.off()
airdata_ntx_all_all <- subset(airdata_ntx_all,type=="all")
wilcox.test(delta_p_pt~gn,data=airdata_ntx_all_all)
airdata_ntx_all_art <- subset(airdata_ntx_all,type=="a")
wilcox.test(delta_p_pt~gn,data=airdata_ntx_all_art)
airdata_ntx_all_cap <- subset(airdata_ntx_all,type=="c")
wilcox.test(delta_p_pt~gn,data=airdata_ntx_all_cap)
airdata_ntx_all_ven <- subset(airdata_ntx_all,type=="v")
wilcox.test(delta_p_pt~gn,data=airdata_ntx_all_ven)


## CAA
# airdata_matched_ntx_caa contains rows with diving vessels that have CAA
# airdata_matched_ntx_pial contains rows with pial vessels that have CAA
# airdata_matched_ntx_plaque contain rows of vessels from imaging window that has plaque
data_matched_ntx <- subset(data_matched,trt=="ntx")
data_matched_ntx_caa <- data_matched_ntx[complete.cases(data_matched_ntx[,"vessel_caa"]),]
data_matched_ntx_caa$v_caa_density <- data_matched_ntx_caa$vessel_caa/data_matched_ntx_caa$vessel_area
data_matched_ntx_caa$p_caa_density <- data_matched_ntx_caa$pial_caa/data_matched_ntx_caa$pial_area
airdata_matched_ntx_caa <- subset(data_matched_ntx_caa,state=="air")
cdata_matched_ntx_caa <- subset(data_matched_ntx_caa,state=="co2")
delta_pt <- (cdata_matched_ntx_caa$pt - airdata_matched_ntx_caa$pt)/airdata_matched_ntx_caa$pt
delta_slope = (cdata_matched_ntx_caa$slope - airdata_matched_ntx_caa$slope)*100/airdata_matched_ntx_caa$slope
delta_AUC = (cdata_matched_ntx_caa$AUC - airdata_matched_ntx_caa$AUC)*100/airdata_matched_ntx_caa$AUC
airdata_matched_ntx_caa$delta_pt <- delta_pt
airdata_matched_ntx_caa$delta_slope <- delta_slope
airdata_matched_ntx_caa$delta_AUC <- delta_AUC
airdata_matched_ntx_caa <- subset(airdata_matched_ntx_caa,v_caa_density<0.72)
cdata_matched_ntx_caa$delta_pt <- delta_pt
cdata_matched_ntx_caa$delta_slope <- delta_slope
cdata_matched_ntx_caa$delta_AUC <- delta_AUC
cdata_matched_ntx_caa <- subset(cdata_matched_ntx_caa,v_caa_density<0.72)
data_matched_ntx_caa <- rbind(airdata_matched_ntx_caa,cdata_matched_ntx_caa)

data_matched_ntx_tg <- subset(data_matched_ntx,gn=="tg"&type!="c")
data_matched_ntx_tg$v_caa_density <- data_matched_ntx_tg$vessel_caa/data_matched_ntx_tg$vessel_area
data_matched_ntx_tg$P_caa_density <- data_matched_ntx_tg$pial_caa/data_matched_ntx_tg$pial_area
airdata_matched_ntx_tg <- subset(data_matched_ntx_tg,state=="air")
cdata_matched_ntx_tg <- subset(data_matched_ntx_tg,state=="co2")
delta_pt <- airdata_matched_ntx_tg$pt - cdata_matched_ntx_tg$pt
delta_p_pt <- (airdata_matched_ntx_tg$pt - cdata_matched_ntx_tg$pt)/airdata_matched_ntx_tg$pt*100
delta_slope = (cdata_matched_ntx_tg$slope - airdata_matched_ntx_tg$slope)*100/airdata_matched_ntx_tg$slope
delta_AUC = (cdata_matched_ntx_tg$AUC - airdata_matched_ntx_tg$AUC)*100/airdata_matched_ntx_tg$AUC
airdata_matched_ntx_tg$delta_pt <- delta_pt
airdata_matched_ntx_tg$delta_p_pt <- delta_p_pt
airdata_matched_ntx_tg$delta_slope <- delta_slope
airdata_matched_ntx_tg$delta_AUC <- delta_AUC
airdata_matched_ntx_tg$v_caa_density[is.na(airdata_matched_ntx_tg$v_caa_density)] <- 0
#v_caa_density VS change in TTP PLOT
# airdata_matched_ntx_tg <-subset(airdata_matched_ntx_tg,type=="a")
# caa_ac <-data.frame(airdata_matched_ntx_tg$v_caa_density,airdata_matched_ntx_tg$delta_pt)
# colnames(caa_ac) <- c("v_caa_density","delta_pt")
bb = ggplot(data=airdata_matched_ntx_tg,aes(x=v_caa_density,y=delta_p_pt,shape=type))+geom_point(size=6)+
  ggtitle("")+xlim(0,0.40)+
  scale_shape_discrete(guide = FALSE)+theme_gray()+
  scale_shape_manual("Vessel Type",values=c(16,2),labels = c("Arterioles","Venules")) + 
#   labs(x="CAA Density",
#        y=expression(atop("Vascular Reactivity", paste("[",TTP[air]," - ",TTP[CO[2]],", s]"))))+
#   labs(x="CAA Density",
#        y=expression(atop("Vascular Reactivity", paste("[(",TTP[air]," - ",TTP[CO[2]],")/",TTP[air], ", %]"))))+
    labs(x="CAA Load",y="Vascular Reactivity [%]")+
  theme(axis.text=element_text(size=30),
        text=element_text(size=30),
        legend.position = "bottom")
#         legend.position = c(0.78, 1), 
#         legend.justification = c(0, 1))
bb
dev.copy(jpeg,"vascular_amyloid/V_amyloid_VS_Vascular_reactivity-ART_ONLY-NO_FIT.jpeg", width=800, height=800)
dev.off()
# vertical-distance linear regression
sum = lme(delta_p_pt ~ v_caa_density, random = ~1|id, na.action = na.omit, data=airdata_matched_ntx_tg)
summ = summary(lme(delta_p_pt ~ v_caa_density, random = ~1|id, na.action = na.omit, data=airdata_matched_ntx_tg))
ftd <- fitted(sum)
R2<- 1-sum((airdata_matched_ntx_tg$delta_p_pt-ftd)^2)/sum((airdata_matched_ntx_tg$delta_p_pt-mean(airdata_matched_ntx_tg$delta_p_pt))^2)
R2
cu_slope <- summ$coeff$fixed[2]+sqrt(diag(summ$var))[2]
cl_slope <- summ$coeff$fixed[2]-sqrt(diag(summ$var))[2]
cu_int <- summ$coeff$fixed[1]-sqrt(diag(summ$var))[1]
cl_int <- summ$coeff$fixed[1]+sqrt(diag(summ$var))[1]
airdata_matched_ntx_tg$cu <- cu_slope*airdata_matched_ntx_tg$v_caa_density+cu_int
airdata_matched_ntx_tg$cl <- cl_slope*airdata_matched_ntx_tg$v_caa_density+cl_int
min_x <- min(airdata_matched_ntx_tg$v_caa_density)
max_x <- max(airdata_matched_ntx_tg$v_caa_density)
max_y <- summ$coeff$fixed[1] + summ$coeff$fixed[2]*max_x
lb1 <- paste("R^2 == ", round(R2,2))
# orthogonal-distance linear regression
caa_pca <- deming(airdata_matched_ntx_tg$delta_p_pt~airdata_matched_ntx_tg$v_caa_density)
caa_ci_int <- (caa_pca$ci[2,2]-caa_pca$ci[2,1])/2               # ci interval for slope
caa_slp_se <- sqrt(diag(caa_pca$variance))[2]                   # se for slope
caa_slp_int_se <- sqrt(diag(caa_pca$variance))[1]            # se for original slope's intercept
caa_pca_t <- (caa_pca$coeff[2])/(caa_slp_se)
caa_pca_p <- 2*pt(-abs(caa_pca_t),df=length(airdata_matched_ntx_tg$delta_p_pt)-1)
caa_pca <- deming(airdata_matched_ntx_tg$delta_pt~airdata_matched_ntx_tg$v_caa_density)
caa_ci_int <- (caa_pca$ci[2,2]-caa_pca$ci[2,1])/2               # ci interval for slope
caa_slp_se <- sqrt(diag(caa_pca$variance))[2]                   # se for slope
caa_slp_int_se <- sqrt(diag(caa_pca$variance))[1]            # se for original slope's intercept
caa_pca_t <- (caa_pca$coeff[2])/(caa_slp_se)
caa_pca_p <- 2*pt(-abs(caa_pca_t),df=length(airdata_matched_ntx_tg$delta_pt)-1)
a_load <- bb+geom_segment(aes(x = min_x, y = summ$coeff$fixed[1], xend = max_x, yend = max_y[1]),
                      size=1.5,color="black")
#          stat_smooth(data=subset(airdata_matched_ntx_tg,v_caa_density<1),
#                       method = f, se = TRUE, colour = 'red', size=1.5, formula=y~x)
#         annotate("text",x=0.65, y=-0.05, size=12, label= lb1, parse=TRUE)
a_load
dev.copy(jpeg,"vascular_amyloid/adjusted_V_amyloid_VS_Vascular_reactivity-perc_change-orginary_reg.jpeg", width=800, height=800)
dev.off()

blankPlot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()
top_row <- plot_grid(a_load,blankPlot,perc_ttp, labels=c("A","","B"),label_size=25,ncol=3, rel_widths = c(1,0.05,1))
top_row
dev.copy(jpeg,"delta_ttp/a_load&delta_perc_ttp.jpeg", width=1600, height=800)
dev.off()



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

# all_regression wt ntx     NOT CENTRED
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
colnames(wt_ntx_ac) <- c("air","co2")
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
wt_ntx_pca <- deming(wt_ntx_ac[,2]~wt_ntx_ac[,1])
wt_ntx_flow_change <- 1/wt_ntx_pca$coefficients[2]
wt_ntx_ci_int <- (wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1])/2
wt_ntx_flow_change_ci <- wt_ntx_ci_int*(wt_ntx_flow_change)^2
wt_ntx_slp_se <- sqrt(diag(wt_ntx_pca$variance))[2]                     # se for slope
wt_ntx_flow_se <- wt_ntx_slp_se/wt_ntx_pca$coefficients[2]^2            # se for flow change
wt_ntx_slp_int_se <- sqrt(diag(wt_ntx_pca$variance))[1]                 # se for original slope's intercept
wt_ntx_flow_int_se <- wt_ntx_slp_int_se/wt_ntx_pca$coefficients[1]^2    # se for flow change
line_cu <- list(slope=wt_ntx_pca$ci[2,2],intercept=wt_ntx_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_pca$ci[2,1],intercept=wt_ntx_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_ac$se_u <- (wt_ntx_pca$coefficients[2]+wt_ntx_slp_se)*wt_ntx_ac$air+wt_ntx_pca$coeff[1]-wt_ntx_slp_int_se
wt_ntx_ac$se_l <- (wt_ntx_pca$coefficients[2]-wt_ntx_slp_se)*wt_ntx_ac$air+wt_ntx_pca$coeff[1]+wt_ntx_slp_int_se
wt_ntx_ac$cu <- wt_ntx_pca$ci[2,2]*wt_ntx_ac$air+wt_ntx_pca$ci[1,1]
wt_ntx_ac$cl <- wt_ntx_pca$ci[2,1]*wt_ntx_ac$air+wt_ntx_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_pca$ci[1,2]-wt_ntx_pca$ci[1,1])/(wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1]),
                             y=(wt_ntx_pca$ci[2,1]*wt_ntx_pca$ci[1,1]-wt_ntx_pca$ci[2,2]*wt_ntx_pca$ci[1,2])/(wt_ntx_pca$ci[2,1]-wt_ntx_pca$ci[2,2])))
wt_ntx_plot <- ggplotRegressionPCA(wt_ntx_ac,5.5,12)
wt_ntx_plot
aa <- wt_ntx_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("nTg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_nTG_nTX_ALL.jpeg")
dev.off()
wt_ntx_flow_change_tab = paste(round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_se,digits=3))
# all_regression tg ntx     NOT CENTRED
tg_ntx_ac <-data.frame(airdata_tg_ntx_mtchd$pt,cdata_tg_ntx_mtchd$pt)
colnames(tg_ntx_ac) <- c("air","co2")
tg_ntx_pca <- deming(tg_ntx_ac[,2]~tg_ntx_ac[,1])
tg_ntx_flow_change <- 1/tg_ntx_pca$coefficients[2]
tg_ntx_ci_int <- (tg_ntx_pca$ci[2,2]-tg_ntx_pca$ci[2,1])/2              # ci interval for slope
tg_ntx_flow_change_ci <- tg_ntx_ci_int*(tg_ntx_flow_change)^2           # ci interval for flow change
tg_ntx_slp_se <- sqrt(diag(tg_ntx_pca$variance))[2]                     # se for slope
tg_ntx_flow_se <- tg_ntx_slp_se/tg_ntx_pca$coefficients[2]^2            # se for flow change
tg_ntx_slp_int_se <- sqrt(diag(tg_ntx_pca$variance))[1]                 # se for original slope's intercept
tg_ntx_flow_int_se <- tg_ntx_slp_int_se/tg_ntx_pca$coefficients[1]^2    # se for flow change
tg_ntx_ac$se_u <- (tg_ntx_pca$coefficients[2]+tg_ntx_slp_se)*tg_ntx_ac$air+tg_ntx_pca$coeff[1]-tg_ntx_slp_int_se
tg_ntx_ac$se_l <- (tg_ntx_pca$coefficients[2]-tg_ntx_slp_se)*tg_ntx_ac$air+tg_ntx_pca$coeff[1]+tg_ntx_slp_int_se
tg_ntx_ac$cu <- tg_ntx_pca$ci[2,2]*tg_ntx_ac$air+tg_ntx_pca$ci[1,1]
tg_ntx_ac$cl <- tg_ntx_pca$ci[2,1]*tg_ntx_ac$air+tg_ntx_pca$ci[1,2]
tg_ntx_plot <- ggplotRegressionPCA(tg_ntx_ac,5,12)
aa <- tg_ntx_plot+xlim(4,12.5)+ylim(4,12.5)+
  labs(y=expression(paste("TTP during ",CO[2]," [s]")))+labs(x="TTP during air [s]")+ggtitle("Tg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
dev.copy(jpeg,"FC_TG_nTX_ALL.jpeg")
dev.off()
tg_ntx_flow_change_tab = paste(round(tg_ntx_flow_change,digits=3)," +/- ",round(tg_ntx_flow_se,digits=3))
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

# all_regression wt ntx     CENTRED
wt_ntx_ac <-data.frame(airdata_wt_ntx_mtchd$pt,cdata_wt_ntx_mtchd$pt)
colnames(wt_ntx_ac) <- c("air","co2")
wt_ntx_ac$TP_air_centered <- wt_ntx_ac$air - mean(wt_ntx_ac$air,na.rm=TRUE)
wt_ntx_ac$TP_co2_centered <- wt_ntx_ac$co2 - mean(wt_ntx_ac$co2,na.rm=TRUE)
ggplotRegressionPCA <- function (df,x_left,x_right) {
  require(ggplot2)
  ggplot(df,aes_string(names(df)[3],names(df)[4]))+geom_point()+
    geom_line(aes(x=TP_air_centered,y=cu),colour="darkgrey",alpha="0.5")+
    geom_line(aes(x=TP_air_centered,y=cl),colour="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cl,ymax=cu),fill="darkgrey",alpha="0.5")+
    geom_ribbon(aes(ymin=cu,ymax=cl),fill="darkgrey",alpha="0.5")+
    stat_smooth(data=subset(df,TP_co2_centered>x_left&TP_co2_centered<x_right),
                method = f, se = TRUE, colour = 'red', size=1.5, formula=y~x)
}
wt_ntx_pca <- deming(wt_ntx_ac[,4]~wt_ntx_ac[,3])
wt_ntx_flow_change <- 1/wt_ntx_pca$coefficients[2]
wt_ntx_ci_int <- (wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1])/2
wt_ntx_flow_change_ci <- wt_ntx_ci_int*(wt_ntx_flow_change)^2
wt_ntx_slp_se <- sqrt(diag(wt_ntx_pca$variance))[2]                     # se for slope
wt_ntx_flow_se <- wt_ntx_slp_se/wt_ntx_pca$coefficients[2]^2            # se for flow change
wt_ntx_slp_int_se <- sqrt(diag(wt_ntx_pca$variance))[1]                 # se for original slope's intercept
wt_ntx_flow_int_se <- wt_ntx_slp_int_se/wt_ntx_pca$coefficients[1]^2    # se for flow change
line_cu <- list(slope=wt_ntx_pca$ci[2,2],intercept=wt_ntx_pca$ci[1,1])
line_cl <- list(slope=wt_ntx_pca$ci[2,1],intercept=wt_ntx_pca$ci[1,2])
lines <- as.data.frame(seq(5.65,11.65,by=0.1))
colnames(lines) <- "x"
wt_ntx_ac$se_u <- (wt_ntx_pca$coefficients[2]+wt_ntx_slp_se)*wt_ntx_ac$TP_air_centered+wt_ntx_pca$coeff[1]-wt_ntx_slp_int_se
wt_ntx_ac$se_l <- (wt_ntx_pca$coefficients[2]-wt_ntx_slp_se)*wt_ntx_ac$TP_air_centered+wt_ntx_pca$coeff[1]+wt_ntx_slp_int_se
wt_ntx_ac$cu <- wt_ntx_pca$ci[2,2]*wt_ntx_ac$TP_air_centered+wt_ntx_pca$ci[1,1]
wt_ntx_ac$cl <- wt_ntx_pca$ci[2,1]*wt_ntx_ac$TP_air_centered+wt_ntx_pca$ci[1,2]
intersect <- as.data.frame(c(x=(wt_ntx_pca$ci[1,2]-wt_ntx_pca$ci[1,1])/(wt_ntx_pca$ci[2,2]-wt_ntx_pca$ci[2,1]),
                             y=(wt_ntx_pca$ci[2,1]*wt_ntx_pca$ci[1,1]-wt_ntx_pca$ci[2,2]*wt_ntx_pca$ci[1,2])/(wt_ntx_pca$ci[2,1]-wt_ntx_pca$ci[2,2])))
wt_ntx_plot <- ggplotRegressionPCA(wt_ntx_ac,-3,3.5)
wt_ntx_plot
aa <- wt_ntx_plot+xlim(-3,3.5)+ylim(-3,3.5)+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("nTg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
# dev.copy(jpeg,"FC_nTG_nTX_ALL.jpeg")
# dev.off()
wt_ntx_flow_change_tab = paste(round(wt_ntx_flow_change,digits=3)," +/- ",round(wt_ntx_flow_se,digits=3))
# all_regression tg ntx     CENTRED
tg_ntx_ac <-data.frame(airdata_tg_ntx_mtchd$pt,cdata_tg_ntx_mtchd$pt)
colnames(tg_ntx_ac) <- c("air","co2")
tg_ntx_ac$TP_air_centered <- tg_ntx_ac$air - mean(tg_ntx_ac$air,na.rm=TRUE)
tg_ntx_ac$TP_co2_centered <- tg_ntx_ac$co2 - mean(tg_ntx_ac$co2,na.rm=TRUE)
tg_ntx_pca <- deming(tg_ntx_ac[,4]~tg_ntx_ac[,3])
tg_ntx_flow_change <- 1/tg_ntx_pca$coefficients[2]
tg_ntx_ci_int <- (tg_ntx_pca$ci[2,2]-tg_ntx_pca$ci[2,1])/2              # ci interval for slope
tg_ntx_flow_change_ci <- tg_ntx_ci_int*(tg_ntx_flow_change)^2           # ci interval for flow change
tg_ntx_slp_se <- sqrt(diag(tg_ntx_pca$variance))[2]                     # se for slope
tg_ntx_flow_se <- tg_ntx_slp_se/tg_ntx_pca$coefficients[2]^2            # se for flow change
tg_ntx_slp_int_se <- sqrt(diag(tg_ntx_pca$variance))[1]                 # se for original slope's intercept
tg_ntx_flow_int_se <- tg_ntx_slp_int_se/tg_ntx_pca$coefficients[1]^2    # se for flow change
tg_ntx_ac$se_u <- (tg_ntx_pca$coefficients[2]+tg_ntx_slp_se)*tg_ntx_ac$TP_air_centered+tg_ntx_pca$coeff[1]-tg_ntx_slp_int_se
tg_ntx_ac$se_l <- (tg_ntx_pca$coefficients[2]-tg_ntx_slp_se)*tg_ntx_ac$TP_air_centered+tg_ntx_pca$coeff[1]+tg_ntx_slp_int_se
tg_ntx_ac$cu <- tg_ntx_pca$ci[2,2]*tg_ntx_ac$TP_air_centered+tg_ntx_pca$ci[1,1]
tg_ntx_ac$cl <- tg_ntx_pca$ci[2,1]*tg_ntx_ac$TP_air_centered+tg_ntx_pca$ci[1,2]
tg_ntx_plot <- ggplotRegressionPCA(tg_ntx_ac,-2.5,4.5)
aa <- tg_ntx_plot+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+ggtitle("Tg - ALL")+
  geom_abline(slope=1,intercept=0,colour="black",linetype="dashed")+
  theme(axis.text=element_text(size=20),text=element_text(size=20))+theme_gray()+
  theme(axis.text=element_text(size=30,colour="white"),
        text=element_text(size=30,color="white"),
        legend.text=element_text(size=30,color="white"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
aa
# dev.copy(jpeg,"FC_TG_nTX_ALL.jpeg")
# dev.off()
tg_ntx_flow_change_tab = paste(round(tg_ntx_flow_change,digits=3)," +/- ",round(tg_ntx_flow_se,digits=3))
# hypothesis testing
#Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859???866.
z_ntx_flow_change = (tg_ntx_pca$coefficient[2]-wt_ntx_pca$coefficient[2])/
  (sqrt(tg_ntx_slp_se^2+wt_ntx_slp_se^2))
p_ntx_flow_change=2*pnorm(-abs(z_ntx_flow_change))

# Combined plot of all vessels of nTg and all vessels of Tg     CENTERED
ggplotRegressionPCAcombined <- function (df1,df2,colorVar1,colorVar2,x_left1,x_right1,x_left2,x_right2) {
  require(ggplot2)
  df3 <- data.frame(gn=c(rep("Tg",nrow(df1)), rep("nTg",nrow(df2))))
  df3 <- cbind(df3,rbind(df1,df2))
  df3 <- cbind(df3,data.frame(matrix(ncol=6,nrow=length(df3[,1]))))
  left_lim <- floor(min(df3$TP_air_centered,df3$TP_co2_centered)/0.5)*0.5
  right_lim <- ceiling(max(df3$TP_air_centered,df3$TP_co2_centered)/0.5)*0.5
  lim <- max(abs(left_lim),right_lim)
  colnames(df3)[(ncol(df3)-5):ncol(df3)] <- c("x","ydiag","ydiagneg","yhoriz","yvert1","yvert2")
  df3[,(ncol(df3)-5)] <- seq(-5,5,length=length(df3[,1]))
  df3[,(ncol(df3)-4)] <- seq(-5,5,length=length(df3[,1]))
  df3[,(ncol(df3)-3)] <- seq(5,-5,length=length(df3[,1]))
  df3[,(ncol(df3)-2)] <- rep(0,length(df3[,1]))
  # below: couldn't figure out how to creat y values for vertical line while following x values, so just created 
  # infinitely small and large values.
  df3[,(ncol(df3)-1)] <- c(rep(-9999999,floor(length(df3[,1])/2)),rep(9999999,ceiling(length(df3[,1])/2)))
  df3[,ncol(df3)] <- c(rep(9999999,floor(length(df3[,1])/2)),rep(-9999999,ceiling(length(df3[,1])/2)))
 
  ggplot(df3,aes_string(names(df3)[4],names(df3)[5]))+geom_point(aes(color=gn),size=2,alpha=0.4)+
    geom_line(data=df3,aes(x=x,y=ydiag),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=df3,aes(x=x,y=ydiagneg),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=df3,aes(x=x,y=yhoriz),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=df3,aes(x=x,y=yvert1),colour="black",linetype="dashed",size=0.8,alpha=0.3)+
    geom_line(data=df3,aes(x=TP_air_centered,y=cu,group=gn),colour="darkgrey",alpha="0.5")+
    geom_line(data=df3,aes(x=TP_air_centered,y=cl,group=gn),colour="darkgrey",alpha="0.5")+
    geom_ribbon(data=df3,aes(x=x,ymin=ydiagneg,ymax=ydiag),fill="olivedrab4",alpha=0.3)+
    geom_ribbon(data=df3,aes(x=x,ymin=ydiag,ymax=yvert1),fill="darkorange3",alpha=0.3)+
    geom_ribbon(data=df3,aes(x=x,ymin=ydiagneg,ymax=yvert2),fill="darkorange3",alpha=0.3)+
    geom_ribbon(data=df3,aes(ymin=cl,ymax=cu,group=gn),fill="darkgrey",alpha=0.6)+
    scale_color_manual(name="Genotype",values=c(colorVar2,colorVar1),
                       labels=c("nTg","Tg"))+
    #     geom_point(data=df2,color=colorVar2,size=2,alpha=0.4)+
    #     geom_line(data=df2,aes(x=air,y=cu),colour="darkgrey",alpha="0.5")+
    #     geom_line(data=df2,aes(x=air,y=cl),colour="darkgrey",alpha="0.5")+
    #     geom_ribbon(data=df2,aes(ymin=cl,ymax=cu),fill="darkgrey",alpha=0.4)+
    #     geom_ribbon(data=df2,aes(ymin=cu,ymax=cl),fill="darkgrey",alpha=0.4)+
    stat_smooth(data=subset(df1,TP_air_centered>x_left1&TP_air_centered<x_right1),
                method = f, se = TRUE, colour = colorVar1, size=1.5, formula=y~x)+
    stat_smooth(data=subset(df2,TP_air_centered>x_left2&TP_air_centered<x_right2),
                method = f, se = TRUE, colour = colorVar2, size=1.5, formula=y~x)+
    coord_cartesian(xlim=c(-lim,lim),ylim=c(-lim,lim))

#     geom_ribbon(data=df3,aes(ymin=yhoriz,ymax=ydiag),fill="olivedrab4",alpha=0.2)+
#     geom_ribbon(data=subset(df3,x<0&x>-3),aes(ymin=ydiag,ymax=ydiagneg),fill="olivedrab4",alpha=0.2)+

}
combined_plot <- ggplotRegressionPCAcombined(tg_ntx_ac,wt_ntx_ac,'darkorange3','olivedrab4',-3,4,-3,4.5)
combined_plot

comb <- combined_plot+
  labs(y=expression(paste(TTP[CO[2]]," centered [s]")))+labs(x=expression(paste(TTP[air]," centered [s]")))+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+theme_gray()+
  theme(axis.text=element_text(size=30),text=element_text(size=30),
        legend.position="none")
comb
dev.copy(jpeg,"flow_change/TTPair_VS_TTPco2_bg.jpeg", width=800, height=800)
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
# dev.copy(jpeg,"FC_nTG_nTX_Arterioles.jpeg")
# dev.off()
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
# dev.copy(jpeg,"FC_TG_nTX_Arterioles.jpeg")
# dev.off()
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
# dev.copy(jpeg,"FC_nTG_nTX_Capillaries.jpeg")
# dev.off()
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
# dev.copy(jpeg,"FC_TG_nTX_Capillaries.jpeg")
# dev.off()
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
# dev.copy(jpeg,"FC_nTG_nTX_Venules.jpeg")
# dev.off()
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
# dev.copy(jpeg,"FC_TG_nTX_Venules.jpeg")
# dev.off()
tg_ntx_ven_flow_change_tab = paste(round(tg_ntx_ven_flow_change,digits=3)," +/- ",round(tg_ntx_ven_flow_se,digits=3))
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

# % change in flow -> (RFC-1)*100 %
rfc.df <- data.frame(matrix(ncol=4,nrow=8))
colnames(rfc.df) <- c("type","gn","rfc","se")
rfc.df[,3]<-c(wt_ntx_flow_change,tg_ntx_flow_change,wt_ntx_art_flow_change,tg_ntx_art_flow_change,
              wt_ntx_cap_flow_change,tg_ntx_cap_flow_change,wt_ntx_ven_flow_change,tg_ntx_ven_flow_change)
rfc.df[,4]<-c(wt_ntx_flow_se,tg_ntx_flow_se,wt_ntx_art_flow_se,tg_ntx_art_flow_se,
              wt_ntx_cap_flow_se,tg_ntx_cap_flow_se,wt_ntx_ven_flow_se,tg_ntx_ven_flow_se)*100
# rfc.df[,4]<-(rfc.df[,4]/(rfc.df[,3]-1))*100 
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
  geom_errorbar(data=subset(rfc.df,rfc<0),aes(fill=gn,ymin=rfc, ymax=rfc-se),
                size=.7,
                width=.4,
                position=position_dodge(0.9))+
  geom_bar(aes(fill=gn, width=width),position=position_dodge(), stat="identity",
             color="black",
             size=1)+ylim(-25,145)+xlab("")+ylab(expression(paste(Delta,"Flow [%]")))+
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
                        rfc = c(110,124,137))
data1 <- data.frame(x = c(2,2,4,4), y = c(100,110,110,100),
                    type=NA, delta_pt=NA)
data2 <- data.frame(x = c(6.5,6.5,7.5,7.5), y = c(114,124,124,114),
                    type=NA, delta_pt=NA)
data3 <- data.frame(x = c(9.5,9.5,10.5,10.5), y = c(127,137,137,127),
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



## COV TTP
data_ntx = subset(data_matched,trt=="ntx")
# data_count <- data.frame(matrix(ncol=4))
# colnames(data_count) <- c("id","loc","gn","count")
# count <- 1
# for (i in unique(data_ntx$id)) {
#   for (k in unique(data_ntx[data_ntx$id==i,]$loc)) {
#     num_dt <- length(data_ntx[data_ntx$id==i&data_ntx$loc==k,]$pt)
#     if (num_dt < 5) {
#     data_count[count,1] = as.numeric(as.character(data_ntx[data_ntx$id==i,]$id[1]))
#     data_count[count,2] = as.numeric(as.character(data_ntx[data_ntx$loc==k,]$loc[1]))
#     data_count[count,3] = as.character(data_ntx[data_ntx$id==i,]$gn[1])
#     data_count[count,4] = num_dt
#     count <- count+1}
#   }
# }
# for (i in unique(data_count$id)) {
#   location <- data_count[data_count$id==i,]$loc
#   data_ntx<-data_ntx[!(data_ntx$id==i & data_ntx$loc==location),]
# }
# # check that the new 'data_count' is empty
# new_data_count <- data.frame(matrix(ncol=4))
# colnames(new_data_count) <- c("id","loc","gn","count")
# count <- 1
# for (i in unique(data_ntx$id)) {
#   for (k in unique(data_ntx[data_ntx$id==i,]$loc)) {
#     num_dt <- length(data_ntx[data_ntx$id==i&data_ntx$loc==k,]$pt)
#     if (num_dt < 5) {
#       new_data_count[count,1] = as.numeric(as.character(data_ntx[data_ntx$id==i,]$id[1]))
#       new_data_count[count,2] = as.numeric(as.character(data_ntx[data_ntx$loc==k,]$loc[1]))
#       new_data_count[count,3] = as.character(data_ntx[data_ntx$id==i,]$gn[1])
#       new_data_count[count,4] = num_dt
#       count <- count+1}
#   }
# }
# new_data_count
airdata_ntx_wt <- subset(data_ntx,state=="air"&gn=="wt")
cdata_ntx_wt <- subset(data_ntx,state=="co2"&gn=="wt")
airdata_ntx_tg <- subset(data_ntx,state=="air"&gn=="tg")
cdata_ntx_tg <- subset(data_ntx,state=="co2"&gn=="tg")
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
colnames(cov_all) <- c("art_air","art_co2","cap_air","cap_co2","ven_air","ven_co2")
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="a",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="c",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_ntx_wt[airdata_ntx_wt$type=="v",]$pt,cdata_ntx_wt[cdata_ntx_wt$type=="v",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="a",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="a",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="c",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="c",]$pt,alternative=c("greater"))
var.test(airdata_ntx_tg[airdata_ntx_tg$type=="v",]$pt,cdata_ntx_tg[cdata_ntx_tg$type=="v",]$pt,alternative=c("greater"))

sd_id_air_wt <- unclass(by(airdata_ntx_wt[,"pt"],airdata_ntx_wt$id,sd))
sd_id_co2_wt <- unclass(by(cdata_ntx_wt[,"pt"],cdata_ntx_wt$id,sd))
sd_id_air_tg <- unclass(by(airdata_ntx_tg[,"pt"],airdata_ntx_tg$id,sd))
sd_id_co2_tg <- unclass(by(cdata_ntx_tg[,"pt"],cdata_ntx_tg$id,sd))
df_air_wt <- as.data.frame(unclass(sd_id_air_wt))
df_air_wt <- data.frame(df_air_wt[complete.cases(df_air_wt),])
df_co2_wt <- as.data.frame(unclass(sd_id_co2_wt))
df_co2_wt <- data.frame(df_co2_wt[complete.cases(df_co2_wt),])
df_air_tg <- as.data.frame(unclass(sd_id_air_tg))
df_air_tg <- data.frame(df_air_tg[complete.cases(df_air_tg),])
df_co2_tg <- as.data.frame(unclass(sd_id_co2_tg))
df_co2_tg <- data.frame(df_co2_tg[complete.cases(df_co2_tg),])
df_ratio_wt <- df_air_wt/df_co2_wt
df_ratio_tg <- df_air_tg/df_co2_tg
df_ratio_gn <- data.frame(id=c(rownames(df_air_wt),rownames(df_air_tg)),
                          gn= c(rep("wt",nrow(df_ratio_wt)),rep("tg",nrow(df_ratio_tg))),
                          sd_ratio = c(df_ratio_wt[,1],df_ratio_tg[,1]))
df_ratio_gn <- df_ratio_gn[order(rev(df_ratio_gn$gn)),]
bb = ggplot(data=df_ratio_gn)+geom_point(aes(x=gn,y=sd_ratio,shape=gn),size=6)+
  geom_text(aes(x=gn,y=sd_ratio,label=id),hjust=2, vjust=0)+
  scale_x_discrete(limits=c("wt","tg"),labels=c("nTg","TgAD"))+
  ggtitle("")+scale_shape_manual(name="Genotype",values=c(19,8),labels=c("TgAD","nTg"))+
  ylab(expression(paste(sigma[TTP[air]],"/",paste(sigma[TTP[CO[2]]]))))+xlab("")+theme_grey()+
  theme(axis.text=element_text(size=40,colour="white"),
        text=element_text(size=40,color="white"),
        legend.text=element_text(size=40,color="white"),
        legend.key= element_rect(fill = "gray42",color = "gray42"),
        plot.background = element_rect(fill = "gray27"),
        panel.background = element_rect(fill = "gray42"),
        legend.background = element_rect(fill = "gray27"),
        legend.key.height=unit(2.5,"line"))
bb
df_ratio_gn_tg <- df_ratio_gn[df_ratio_gn$gn=="tg",]
df_ratio_gn_tg <- df_ratio_gn_tg[df_ratio_gn_tg$sd_ratio<2,]
df_ratio_gn_wt <- df_ratio_gn[df_ratio_gn$gn=="wt",]
df_ratio_gn <- rbind(df_ratio_gn_wt,df_ratio_gn_tg)
summary(lm(sd_ratio ~ gn, na.action = na.omit, df_ratio_gn))
summary(lme(sd_ratio ~ gn, random = ~1|id, na.action = na.omit, data=df_ratio_gn))
t.test(df_ratio_gn[df_ratio_gn$gn=="wt",]$sd_ratio,df_ratio_gn[df_ratio_gn$gn=="tg",]$sd_ratio)
