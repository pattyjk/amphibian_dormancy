##TbCl strandard curve
library(ggplot2)
tbcl<-read.delim('tbcl_std_curves.txt', header=T)

ggplot(tbcl, aes(No_spores, Fluorescens))+
  facet_wrap(~Name)+
  geom_point()+
  theme_bw()+
  ylab("Fluorescens")+
  xlab('# Spores/ml')+
  geom_smooth(method='lm')+
  scale_x_log10()
