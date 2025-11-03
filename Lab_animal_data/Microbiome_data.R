library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
set.seed(1234)

# load in the asv table from the cluster analysis 
asv_table1 <- read.delim("Lab_animal_data/data/asv_table.txt", header=T)
meta1<-read.delim("Lab_animal_data/data/endospore_map.txt", header=T)
meta2<-read.delim("Lab_animal_data/data/endospore_map_run2.txt", header=T)
meta3<-read.delim("Lab_animal_data/data/endospore_map_run3.txt", header=T)
asv_table2 <- read.delim("Lab_animal_data/data/asv_table2.txt", header=T)
asv_table3 <- read.delim("Lab_animal_data/data/asv_table3.txt", header=T)

# function used to Merge all asv tables together 
MyMerge <- function(x, y){
  df <- merge(x, y, by= "X.OTU.ID", all.x= TRUE, all.y= TRUE)
  return(df)
}
asv_table <- Reduce(MyMerge, list(asv_table1, asv_table2, asv_table3))

# fix formatting of asv table
asv_table <- column_to_rownames(asv_table, var ="X.OTU.ID")
asv_table <- replace(asv_table, is.na(asv_table), 0)

# merge all metadata together
meta <- bind_rows(meta1, meta2, meta3) |>
  distinct()

#look at sequencing depth
colSums(asv_table)

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=107)
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
#7.6

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color=Species))+
  geom_point(aes(size=2, shape = Nuc_type))+
  #geom_text()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 10.7%")+
  ylab("PC2- 7.6%") +
  stat_ellipse(aes(group = Species), type = "norm", level = 0.95, linetype = "dashed")


#### OTU Richness ####

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(asv_table), sample=107)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(rrarefy(t(asv_table), sample = 107))`)
t.test(larv.alph$Richness, larv.alph$Nuc_type)
# t = 12.147, df = 49, p-value < 2.2e-16

larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(asv_table), sample=107), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)

t.test(larv.alph$Shannon, larv.alph$Nuc_type)
# t = 21.537, df = 49, p-value < 2.2e-16

#plot richness
ggplot(larv.alph, aes(Species, Richness, fill=Nuc_type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  theme(legend.position = "none")+
  coord_flip()+
  ylab("sOTU Richness")
  #scale_fill_manual(values = c('#f58231', '#4363d8'))

# plot Shannon Diversity
ggplot(larv.alph, aes(Species, Shannon, fill=Nuc_type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")
  #scale_fill_manual(values = c('#f58231', '#4363d8'))

#calculate percent dormant
#split into DNA RNA table
DNA_cols <- grepl("DNA", names(asv_table))
rna_table<-asv_table[, !DNA_cols]
dna_table<-asv_table[, DNA_cols]

#sort to be in same order
dna_table<-dna_table[,sort(names(dna_table))]
rna_table<-rna_table[,sort(names(rna_table))]

#divide tables to get ratios, add one to each value to remove any issues with division
ratio_table<-(rna_table+1)/(dna_table+1)

#calculate percent dormant
above_1 <- colSums(ratio_table >= 1)

# Total number of rows
total_rows <- nrow(ratio_table)

# Calculate the proportion of rows above or equal to 1 for each column
results <- above_1 / total_rows

# Create a new data frame to store the results
lab_dorm <- data.frame(Column = names(results), Proportion = results)

#append metadata for stuff
lab_dorm<-merge(lab_dorm, meta, by.x='Column', by.y='SampleID')

#plot it
ggplot(lab_dorm, aes(Species, Proportion))+
  geom_boxplot()+
  theme_bw()+
  ylab("Percent Active")+
  xlab("")+
  coord_flip()
