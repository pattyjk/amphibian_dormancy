#' ----------------------------------------------
#' 
#' Andrew Miller-Klugman
#' Lab animal Microbiome Analysis 
#' 
#' ----------------------------------------------

library(vegan)
library(ggplot2)

# load in the asv table from the cluster analysis 
asv_table <- read.delim("Lab_animal_data/data/asv_table.txt", row.names=1, header=T)
meta<-read.delim("Lab_animal_data/data/endospore_map.txt", header=T, row.names = (1))

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
#10.6
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#8

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, fill = Nuc_type))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 10.6%")+
  ylab("PC2- 8%")


#### TEST WTF AM I DOUNG ####

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(asv_table), sample=107)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(rrarefy(t(asv_table), sample = 107))`)
t.test(larv.alph$Richness, larv.alph$Type2)
#t = 11.267, df = 29, p-value < 4.122e-12

larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(asv_table), sample=107), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)

t.test(larv.alph$Shannon, larv.alph$Type2)
# t = 15.044, df = 29, p-value < 3.116e-15

#plot richness
ggplot(larv.alph, aes(Nuc_type, Richness, fill=Nuc_type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTU Richness")+
  scale_fill_manual(values = c('#f58231', '#4363d8'))



#### TESTINGGGGGGGGGGG ####
install.packages("BiocManager")
BiocManager::install("phyloseq", force = T)
library(phyloseq)
library(rlan)

physeq <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE), sample_data(meta))
richness_vector <- estimate_richness(physeq, measures = "observed")




