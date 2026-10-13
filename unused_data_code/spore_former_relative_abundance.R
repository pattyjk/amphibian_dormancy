library(stringi)
library(reshape2)
library(ggpubr)
library(dplyr)
library(plyr)
library(ggplot2)

#read data
library(readr)
lab_taxonomy1 <- read.csv("Lab_animal_data/data/lab_animals_level7.csv", header=T, row.names = 1)
lab_taxonomy2 <- read.csv("Lab_animal_data/data/lab_animals_level7(2).csv", header=T, row.names = 1)
lab_taxonomy<-t(lab_taxonomy)
spore_genera<-read.delim("spore_forming_genera.txt", header=F)

#change to relative abundance
lab_taxonomy<-as.data.frame(sweep(lab_taxonomy, 2, colSums(lab_taxonomy), '/'))

#split taxonomy column to select genera
lab_taxonomy$tax<-row.names(lab_taxonomy)
library(tidyr)
lab_taxonomy<-separate(lab_taxonomy, tax, sep='\\.', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
spore_tax<-lab_taxonomy[lab_taxonomy$Genus %in% spore_genera$V1,]

#reshape data for plotting
library(reshape2)
spore_m<-melt(spore_tax)

#read in metadata
meta<-read.csv("Lab_animal_data/data/sample_data_meta.csv", header=T)

#merge metadata to spore data
spore_m<-merge(spore_m, meta, by.x='variable', by.y='SampleID')

#plot all samples
ggplot(spore_m, aes(variable, value, fill=Genus))+
  theme_bw()+
  ylab('Relative Abundance')+
  xlab("")+
  geom_bar(stat='identity')

#calcualte average per species/nucleic acid type
library(dplyr)
spore_sum<-ddply(spore_m, c("Species.y", "Nuc_type", "Genus"), summarize, mean=mean(value))

#plot mean per species
ggplot(spore_sum, aes(Species.y, mean, fill=Genus))+
  theme_bw()+
  ylab('Relative Abundance')+
  facet_wrap(~Nuc_type)+
  xlab("")+
  geom_bar(stat='identity')
#lol this makes a middle finger :)  

##############################################
#Ribbitr data

#read data
library(readr)
lab_taxonomy <- read.csv("Ribbitr_data/data/panama_level7.csv", header=T, row.names = 1)
lab_taxonomy<-t(lab_taxonomy)
spore_genera<-read.delim("spore_forming_genera.txt", header=F)

#change to relative abundance
lab_taxonomy<-as.data.frame(sweep(lab_taxonomy, 2, colSums(lab_taxonomy), '/'))

#split taxonomy column to select genera
lab_taxonomy$tax<-row.names(lab_taxonomy)
library(tidyr)
lab_taxonomy<-separate(lab_taxonomy, tax, sep='\\.', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
spore_tax<-lab_taxonomy[lab_taxonomy$Genus %in% spore_genera$V1,]

#reshape data for plotting
library(reshape2)
spore_m<-melt(spore_tax)

#read in metadata
meta<-read.delim("Ribbitr_data/data/2022_Panama_16S_map.txt", header=T)

#merge metadata to spore data
spore_m<-merge(spore_m, meta, by.x='variable', by.y='SampleID')

#plot all samples
ggplot(spore_m, aes(variable, value, fill=Genus))+
  theme_bw()+
  ylab('Relative Abundance')+
  xlab("")+
  geom_bar(stat='identity')

#calcualte average per species/nucleic acid type
library(dplyr)
spore_sum<-ddply(spore_m, c("Species.y", "Type", "Genus", 'region', 'site'), summarize, mean=mean(value))

#plot mean per species
ggplot(spore_sum, aes(Species.y, mean, fill=Genus))+
  theme_bw()+
  ylab('Relative Abundance')+
  facet_wrap(~Type)+
  xlab("")+
  geom_bar(stat='identity')
#lol this makes a middle finger too :)