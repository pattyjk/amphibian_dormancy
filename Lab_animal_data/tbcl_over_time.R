#TbCl timecourse analysis

library(ggplot2)
library(reshape2)
tbcl_timecourse <- read.delim("~/GitHub/amphibian_dormancy/Lab_animal_data/data/tbcl_timecourse.txt")

#reshape data for plotting
tbcl_timecourse<-tbcl_timecourse[,-4]
tbcl_tc_m<-melt(tbcl_timecourse)

#rename time column (here, its called variable)
tbcl_tc_m$variable<-gsub('X24hr', '24', tbcl_tc_m$variable)
tbcl_tc_m$variable<-gsub('X16hr', '16', tbcl_tc_m$variable)
tbcl_tc_m$variable<-gsub('X4hr', '4', tbcl_tc_m$variable)

#plot data, facet by species
ggplot(tbcl_tc_m, aes(reorder(variable, sort(as.numeric(variable))), value))+
  geom_boxplot()+
  facet_wrap(~Species)+
  theme_bw()+
  ylab("Fluorescens")+
  xlab("Time (hrs)")
