#'-------------------------------------------------------------
#'
#'    Andrew Miller-Klugman
#'    RIBBiTR Microbiome Analysis
#'    
#'--------------------------------------------------------------

library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
set.seed(1234)

# Read in the asv files and metadata
asv_table1 <- read.delim("Ribbitr_data/data/PaBr_asv_table.txt", header =T)
asv_table2 <- read.delim("Ribbitr_data/data/Panama_asv_table.txt", header =T)
meta1 <-  read.delim("Ribbitr_data/data/2022_PaBr_16s_map.txt", header =T)
meta2 <-  read.delim("Ribbitr_data/data/2022_Panama_16s_map.txt", header =T)

# function used to Merge all asv tables together 
MyMerge <- function(x, y){
  df <- merge(x, y, by= "OTU.ID", all.x= TRUE, all.y= TRUE)
  return(df)
}
asv_table <- Reduce(MyMerge, list(asv_table1, asv_table2))

# fix formatting of asv table
asv_table <- column_to_rownames(asv_table, var ="OTU.ID")
asv_table <- replace(asv_table, is.na(asv_table), 0)

# merge all metadata together
meta <- bind_rows(meta1, meta2) |>
  distinct()

#### PCoA plotting ####

#look at sequencing depth
colSums(asv_table)

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=511)
# loss of 3 samples 

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(nut_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#10.7
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#6

# Removes data that does not contain a value
ko.coords_loc <- ko.coords |>
  filter(Location != "")

#plot PCoA by location 
ggplot(ko.coords_loc, aes(MDS1, MDS2, color=Location, shape=Type))+
  geom_point(aes(size=2, shape = Type))+
  #geom_text()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 10.8%")+
  ylab("PC2- 6%") +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95, linetype = "dashed")

ko.coords_species <- ko.coords |>
 filter(Species != "N/A") |>
  filter(Species != "Negative Control") |>
  filter(Species != "Negative control") |>
  filter(Species != "")

#plot PCoA by species
ggplot(ko.coords_species, aes(MDS1, MDS2, color=Species, shape=Type))+
  geom_point(aes(size=2, shape = Type))+
  #geom_text()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 10.7%")+
  ylab("PC2- 6%") +
  stat_ellipse(aes(group = Species), type = "norm", level = 0.95, linetype = "dashed")


#### OTU Richness ####
#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(asv_table), sample=511)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(rrarefy(t(asv_table), sample = 511))`)
t.test(larv.alph$Richness, larv.alph$Type2)
# t = 12.147, df = 49, p-value < 2.2e-16

larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(asv_table), sample=511), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)

t.test(larv.alph$Shannon, larv.alph$Type2)
# t = 40.079, df = 462, p-value < 2.2e-16

# removes N/a and negatiove controls 
larv_alph_spec <- larv.alph |>
  filter(Species != "N/A") |>
  filter(Species != "Negative Control") |>
  filter(Species != "Negative control") |>
  filter(Species != "")


#plot richness
ggplot(larv_alph_spec, aes(Species, Richness, fill=Location))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  facet_wrap(~Type)+
  ylab("sOTU Richness")
#scale_fill_manual(values = c('#f58231', '#4363d8'))

# plot Shannon Diversity
ggplot(larv_alph_spec, aes(Species, Shannon, fill=Location))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  facet_wrap(~Type)+
  ylab("Shannon Diversity")
#scale_fill_manual(values = c('#f58231', '#4363d8'))

#write alphpa diversity to file
write.table(larv_alph_spec, 'ribbitr_alpha_div.txt', sep='\t', quote=F, row.names=F)

####Calculate percent dormant based on 16S ratios, rarefied data
nut_rare2<-as.data.frame(t(nut_rare))

#create RNA/DNA only tables for calculations
library(dplyr)
pan_rna<-select(nut_rare2, contains(".RNA"))
pan_dna<-select(nut_rare2, contains(".DNA"))

#sort columns to be in the same order
pan_dna<-pan_dna[,order(colnames(pan_dna))]
pan_rna<-pan_rna[,order(colnames(pan_rna))]

#figure out which samples don't have matches
#dna_name<-as.data.frame(names(pan_dna))
#rna_name<-as.data.frame(names(pan_rna))
#write.table(dna_name, 'dna.txt', sep='\t', quote=F, row.names = F)
#write.table(rna_name, 'rna.txt', sep='\t', quote=F, row.names = F)

#prune non matching samples
pan_dna<-select(pan_dna, -c(Pan22.DNA, Pan23.DNA, Pan41.DNA, Pan44.DNA, Pan25.DNA,PA97.DNA,	PA81.DNA,	PA43.DNA,	PA205.DNA,	PA16.DNA,	PA194.DNA,	PA147.DNA,	BR5.DNA,	BR41.DNA,	BR13.DNA,	BR14.DNA))
pan_rna<-select(pan_rna, -c(Pan24.RNA, Pan8.RNA, Pan32.RNA, Pan40.RNA, PA44.RNA, Pan39.RNA, Pan8.RNA, PA204.RNA,	BR71.RNA,	BR38.RNA))

#check size
dim(pan_dna)
dim(pan_rna)

#divide tables to get 16S ratio
pan_ratio<-(pan_rna+1)/(pan_dna+1)

#calculate percent dormant
above_1 <- colSums(pan_ratio >= 1)

# Total number of rows
total_rows <- nrow(pan_ratio)

# Calculate the proportion of rows above or equal to 1 for each column
results <- above_1 / total_rows

# Create a new data frame to store the results
pan_dorm <- data.frame(Column = names(results), Proportion = results)

#append metadata for stuff
pan_dorm<-merge(pan_dorm, meta, by.x='Column', by.y='SampleID')

#plot by species/site
ggplot(pan_dorm, aes(Species, Proportion, fill=region))+
  theme_bw()+
  xlab("")+
  geom_boxplot()+
  ylab("Percent active")+
  coord_flip()

