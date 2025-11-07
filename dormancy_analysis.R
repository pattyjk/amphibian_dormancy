##look at correlations between dormancy and metrics
library(ggplot2)
library(ggpubr)

#read in full RIBBiTR dataset
mastah_rib_datah<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/master_ribbitr_meta.txt', header=T)

#filter out controls and RNA
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Type == 'RNA'),]
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Species == 'Negative Control'),]
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Species == 'Negative control'),]

dorm_mean<-ddply(mastah_rib_datah, c('bd_detected', 'Species'), summarize, mean=mean(as.numeric(Per_active)), sd=sd(as.numeric(Per_active)))

ggplot(dorm_mean, aes(bd_detected, mean, color=Species))+
  geom_point()

#based on Bd +/-
ggplot(mastah_rib_datah, aes(bd_detected, as.numeric(Per_active), fill=Species))+
  geom_boxplot()+
  ylab("Percent active bacteria")+
  xlab("Presence of Bd")+
  theme_bw()

#convert Bd load to log10
mastah_rib_datah$log_bd<-log(as.numeric(mastah_rib_datah$bd_mean_its1_copies_per_swab))

#plot against Bd load
ggplot(mastah_rib_datah, aes(as.numeric(bd_mean_its1_copies_per_swab), as.numeric(Per_active)))+
  geom_point()+
  facet_wrap(~Species, scales='free')+
  stat_cor(method='spearman', position='jitter')+
  geom_smooth(method='lm')+
  ylab("Percent active bacteria")+
  xlab("Bd ITS1 copies")+
  theme_bw()
#sig patterns cor CoPa and IsHe