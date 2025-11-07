##loook at correlations between abundance metrics 
library(ggplot2)
library(ggpubr)

#read in full RIBBiTR dataset
mastah_rib_datah<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/master_ribbitr_meta.txt', header=T)

#filter out controls and RNA
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Type == 'RNA' | mastah_rib_datah$Species == 'Negative Control' | mastah_rib_datah$Species == 'Negative control'),]

#make some plots with Spearmans for qPCR and CFUs
ggplot(mastah_rib_datah, aes(as.numeric(S16CopiesAvg), as.numeric(CFU_bact)))+
  geom_point()+
  stat_cor(method = 'spearman')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('CFU/ml')+
  xlab('Numer of 16S rRNA genes')

#good correlations across the board, not great for RaPi and PsCr due to low sample sizes, probably exclude from this

#look at  qPCR and Cell Staining
ggplot(mastah_rib_datah, aes(as.numeric(S16CopiesAvg), as.numeric(Total_cells_ml)))+
  geom_point()+
  stat_cor(method = 'spearman')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('Cells/ml')+
  xlab('Numer of 16S rRNA genes')
#good correlations across the board, not great for RaPi and PsCr due to low sample sizes, probably exclude from this

#look at CFU and cell staining
ggplot(mastah_rib_datah, aes(as.numeric(CFU_bact), as.numeric(Total_cells_ml)))+
  geom_point()+
  stat_cor(method = 'spearman')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('CFU/ml')+
  xlab('Cells/ml')
#not as great correlations across the board, not great for RaPi and PsCr due to low sample sizes, probably exclude from this

##That's cool, does anything else correlate with abundance of cells (only use Total_cells_ml because most accurate and I invested the most time in that)
#By species
ggplot(mastah_rib_datah, aes(Species, as.numeric(Total_cells_ml)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(Total_cells_ml) ~Species, data=mastah_rib_datah))
#Rana's have highest, though RaMu lowest of Rana

#by site
ggplot(mastah_rib_datah, aes(Location, as.numeric(Total_cells_ml)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(Total_cells_ml) ~ Location, data=mastah_rib_datah))
#Penn highest, Br/Pan lowest. So higher latitude=more bacteria

#look at size and sex
ggplot(mastah_rib_datah, aes(as.numeric(snout_vent_length), as.numeric(Total_cells_ml)))+
  geom_point()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('SVL (mm)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#correlates pretty good with size, does this reflect sampling effort or a true effect of size? Looks to be driven by RaCa

ggplot(mastah_rib_datah, aes(as.numeric(body_mass_g), as.numeric(Total_cells_ml)))+
  geom_point()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('Mass (g)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#same as SVL, but again, jsut sampling ease of larger frogs? Looks to be driven by RaCa

ggplot(mastah_rib_datah, aes(as.numeric(total_peptides_ug), as.numeric(Total_cells_ml)))+
  geom_point()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('Total peptides (ug)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#scales with peptide concentrations...weird?

#Bd load?
ggplot(mastah_rib_datah, aes(as.numeric(bd_mean_its1_copies_per_swab), as.numeric(Total_cells_ml)))+
  geom_point()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ylab('Cells/ml')+
  xlab('Bd load (log10)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#no scaling here, so that's cool. But for LiWa by itself, there is a negative/sig correlation

#change based on Bd presence/absence?
ggplot(mastah_rib_datah, aes(bd_detected, as.numeric(Total_cells_ml)))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Species)+
  ylab('Cells/ml')+
  xlab('Bd presence')
#no change in abundance based on +/- Bd on skin

#change based on sex of animal?
ggplot(mastah_rib_datah, aes(sex, as.numeric(Total_cells_ml), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#sex doesn't look important
TukeyHSD(aov(as.numeric(Total_cells_ml) ~ Location + sex, data=mastah_rib_datah))

#life stage?
ggplot(mastah_rib_datah, aes(life_stage, as.numeric(Total_cells_ml), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#Juvs seem to have higher bacterial loads. But that's only Penn frogs. Subadults are a bit higher than adults in Brazil. 
TukeyHSD(aov(as.numeric(Total_cells_ml) ~ Location + life_stage, data=mastah_rib_datah))

##let's look at the abundance of spore-forming bacteria based on CFU counts
#any correlation to total abundance?
ggplot(mastah_rib_datah, aes(as.numeric(CFU_bact), as.numeric(CFU_past)))+
  geom_point()+
  theme_bw()+
  stat_cor(method = 'spearman')+
  ylab('Pasturized CFU/ml')+
  xlab('CFU/ml')+
  facet_wrap(~Species)
#nada here, so maybe patters elsewhere

#By species
ggplot(mastah_rib_datah, aes(Species, as.numeric(CFU_past)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(CFU_past) ~Species, data=mastah_rib_datah))
#Rana's have highest, though RaMu lowest of Rana

#by site
ggplot(mastah_rib_datah, aes(Location, as.numeric(CFU_past)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(CFU_past) ~ Location, data=mastah_rib_datah))
#Penn highest, Br/Pan lowest. So higher latitdue=more spore forming bacteria. That fits

#look at size and sex
ggplot(mastah_rib_datah, aes(as.numeric(snout_vent_length), as.numeric(CFU_past)))+
  geom_point()+
  facet_wrap(~Species)+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('SVL (mm)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Not much, only PsCr but low sample numbers

ggplot(mastah_rib_datah, aes(as.numeric(body_mass_g), as.numeric(CFU_past), color=Species))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Species)+
  ylab('Cells/ml')+
  xlab('Mass (g)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Nada here. 

ggplot(mastah_rib_datah, aes(as.numeric(total_peptides_ug), as.numeric(CFU_past)))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Species, scales='free_x')+
  ylab('Cells/ml')+
  xlab('Total peptides (ug)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Sig Neg correlation for HyPh, but nothing else

ggplot(mastah_rib_datah, aes(as.numeric(bd_mean_its1_copies_per_swab), as.numeric(CFU_past), color=Species))+
  geom_point()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ylab('Cells/ml')+
  facet_wrap(~Species)+
  xlab('Bd load')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Neg sig correlation for RaCa. So more Bd=less spore forming bacteria = less resilience, cool.

ggplot(mastah_rib_datah, aes(bd_detected, as.numeric(CFU_past)))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Species)+
  ylab('Cells/ml')+
  xlab('Bd presence')
#no change in abundance based on +/- Bd on skin

ggplot(mastah_rib_datah, aes(sex, as.numeric(CFU_past), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#Nada, jsut site
TukeyHSD(aov(as.numeric(CFU_past) ~ Location + sex, data=mastah_rib_datah))

ggplot(mastah_rib_datah, aes(life_stage, as.numeric(CFU_past), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#Nada, just site
TukeyHSD(aov(as.numeric(CFU_past) ~ Location + life_stage, data=mastah_rib_datah))
######################################################################################################################################################
##That's culturable bacteria, let's look at spores for the Terbium chloride assay for just total spores, since viability is >99% for all samples
#any correlation to total abundance?
ggplot(mastah_rib_datah, aes(as.numeric(CFU_bact), as.numeric(TotalNoSporesTbClCor)))+
  geom_point()+
  theme_bw()+
  stat_cor(method = 'spearman')+
  ylab('Pasturized CFU/ml')+
  xlab('CFU/ml')+
  facet_wrap(~Species)
#nada here for most. But RaMu (+ cor) and Rapi (-cor) are sig, but RaPi has small samples size. but

#By species
ggplot(mastah_rib_datah, aes(Species, as.numeric(TotalNoSporesTbClCor)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(TotalNoSporesTbClCor) ~Species, data=mastah_rib_datah))
#Rana's have highest, though RaMu highest of all followed by RaPi

#by site
ggplot(mastah_rib_datah, aes(Location, as.numeric(TotalNoSporesTbClCor)))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
TukeyHSD(aov(as.numeric(TotalNoSporesTbClCor) ~ Location, data=mastah_rib_datah))
#California highest

#look at size and sex
ggplot(mastah_rib_datah, aes(as.numeric(snout_vent_length), as.numeric(TotalNoSporesTbClCor)))+
  geom_point()+
  facet_wrap(~Species)+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('SVL (mm)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Nada

ggplot(mastah_rib_datah, aes(as.numeric(body_mass_g), as.numeric(TotalNoSporesTbClCor), color=Species))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Species)+
  ylab('Cells/ml')+
  xlab('Mass (g)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Sig neg for RaCa, similar idea about sampling effort?

ggplot(mastah_rib_datah, aes(as.numeric(total_peptides_ug), as.numeric(TotalNoSporesTbClCor), color=Species))+
  geom_point()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('Total peptides (ug)')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Nada

ggplot(mastah_rib_datah, aes(as.numeric(bd_mean_its1_copies_per_swab), as.numeric(TotalNoSporesTbClCor), color=Species))+
  geom_point()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ylab('Cells/ml')+
  facet_wrap(~Species, scales='free_y')+
  xlab('Bd load')+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')
#Nada

ggplot(mastah_rib_datah, aes(bd_detected, as.numeric(TotalNoSporesTbClCor)))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Species)+
  ylab('Cells/ml')+
  xlab('Bd presence')
#Looks like Bd increases the amount of spores, likely just driven by RaMu
TukeyHSD(aov(as.numeric(TotalNoSporesTbClCor) ~ bd_detected + Species, data=mastah_rib_datah))

ggplot(mastah_rib_datah, aes(sex, as.numeric(TotalNoSporesTbClCor), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#Nada, just site. California leading the way
TukeyHSD(aov(as.numeric(TotalNoSporesTbClCor) ~ Location + sex, data=mastah_rib_datah))

ggplot(mastah_rib_datah, aes(life_stage, as.numeric(TotalNoSporesTbClCor), fill=Location))+
  geom_boxplot()+
  theme_bw()+
  ylab('Cells/ml')+
  xlab('')+
  coord_flip()
#Nada, just site. California leading the way
TukeyHSD(aov(as.numeric(TotalNoSporesTbClCor) ~ Location + life_stage, data=mastah_rib_datah))

