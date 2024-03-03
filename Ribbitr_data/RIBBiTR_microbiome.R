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
  df <- merge(x, y, by= "X.OTU.ID", all.x= TRUE, all.y= TRUE)
  return(df)
}
asv_table <- Reduce(MyMerge, list(asv_table1, asv_table2))

# fix formatting of asv table
asv_table <- column_to_rownames(asv_table, var ="X.OTU.ID")
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
ggplot(ko.coords_loc, aes(MDS1, MDS2, color=Location))+
  geom_point(aes(size=2))+
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
ggplot(ko.coords_species, aes(MDS1, MDS2, color=Species))+
  geom_point(aes(size=2))+
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
  ylab("sOTU Richness")
#scale_fill_manual(values = c('#f58231', '#4363d8'))

# plot Shannon Diversity
ggplot(larv_alph_spec, aes(Species, Shannon, fill=Location))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")
#scale_fill_manual(values = c('#f58231', '#4363d8'))

