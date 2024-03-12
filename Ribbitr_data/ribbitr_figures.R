#----------------------------
#
# Andrew Miller-Klugman
# Ribbitr Data Figures
#
#----------------------------

library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)

# load in the data
cell_counts <- read_xlsx("Ribbitr_data/data/RIBBiTR_CFUs_cell_counts.xlsx", sheet = 1)
TbCl_data <- read_xlsx("Ribbitr_data/data/RIBBiTR_TbCl_data.xlsx")
Panama_cell <- read_xlsx("Ribbitr_data/data/Panama_data.xlsx", sheet = 1)
cell_stains <- read_xlsx("Ribbitr_data/data/RIBBiTR_CFUs_cell_counts.xlsx", sheet = 2)
panama_tbcl <- read_xlsx("Ribbitr_data/data/Panama_data.xlsx", sheet = 3)
panama_counts <-  read_xlsx("Ribbitr_data/data/Panama_data.xlsx", sheet = 2)
ribbitr_cell_count <- read_xlsx("Ribbitr_data/data/RIBBiTR_CFUS.xlsx")

#### Cell Culture Count Data ####

#data cleaning 
names(ribbitr_cell_count)[names(ribbitr_cell_count) == "CFU_bact"] <- "CFU_ml"
names(ribbitr_cell_count)[names(ribbitr_cell_count) == "Percent_spore"] <- "Per_spore" 
panama_counts$Location <- "Panama"

cell_counts_clean <-
  ribbitr_cell_count |>
  select("RIBBiTR_ID", "Species", "CFU_ml", "CFU_past", "Per_spore", "Location")

cell_counts_total <- bind_rows(cell_counts_clean, panama_counts)

# boxplot of bacterial CFUs by species per location
cell_counts_total  |>
  ggplot(aes(Species, CFU_ml, fill = Location)) +
  geom_boxplot()

#boxplot of pasteurized CFUs
cell_counts_total |>
  ggplot( aes(Species, CFU_past, fill = Location)) +
  geom_boxplot()


# Boxplot of %spore forming bacteria 
cell_counts_total |>
  ggplot(aes(Species, (100*(CFU_past/CFU_ml)), fill = Location)) +
  geom_boxplot(outlier.shape = NA, outlier.colour = NA) +
  scale_y_continuous(limits = c(0,10)) +
  ylab("% Culutrable Spore Forming Bacteria")
  

# Percent_spore by month
cell_counts |>
  ggplot(aes(Species...5, Percent_spore, fill = Month_sampled...6)) +
  geom_boxplot(outlier.shape = NA, outlier.colour = NA) +
  scale_y_continuous(limits = c(0,10))

# s16 copies by species
cell_counts |>
  ggplot(aes(Species...5, S16_copies, fill = Location...2)) +
  geom_boxplot()


#### TbCl Data ####

# clean data to make combining possible
panama_tbcl$Location <- "Panama"
panama_tbcl$Per_viable_spore <- (panama_tbcl$Viable_fluor/panama_tbcl$Total_fluor)*100
names(panama_tbcl)[names(panama_tbcl) == "RIBBITR_ID"] <- "RIBBiTR_ID"

TbCl_data_clean <- 
  TbCl_data |>
  select("RIBBiTR_ID", "Species", "Location", "Viable_fluor", "Total_fluor", "Per_viable_spore")
Panama_tbcl_data_clean <- 
  panama_tbcl |>
  select("RIBBiTR_ID", "Species", "Location", "Viable_fluor", "Total_fluor", "Per_viable_spore")

TbCl_data_total <- bind_rows(TbCl_data_clean, Panama_tbcl_data_clean)
#percent viable spore
TbCl_data_total |>
  ggplot(aes(Species, Per_viable_spore, fill = Location)) +
  geom_boxplot()

#### Cell Stain Data ####
# select columns needed from each table to ensure binding together
cell_stain_clean <- cell_stains |>
  select("Per_dormant", "RIBBiTR_ID", "Location", "Species",
         "Sublocation", "Sublocation_2", "Month_sampled",
         "Year_sampled", "Total bacteria_blue", "Total_bacteria_red", 
         "Total_cells_ml", "Active_cells_ml", "Per_active")

Panama_stain <- Panama_cell |>
  select("Per_dormant", "RIBBiTR_ID", "Location", "Species",
         "Sublocation", "Sublocation_2", "Month_sampled",
         "Year_sampled", "Total bacteria_blue", "Total_bacteria_red", 
         "Total_cells_ml", "Active_cells_ml", "Per_active")

# bind both dataframes into a single frame 
cell_stain_counts_total <- bind_rows(cell_stain_clean, Panama_cell)

# convert data to long to allow for graphing on same plot 
cell_stain_count_total_long <-
  pivot_longer(cell_stain_counts_total, cols = c('Total bacteria_blue', 'Total_bacteria_red'), names_to = "stain", values_to = "value")

 cell_stain_count_total_long |>
  ggplot(aes(x = value, y = Species, fill = stain)) +
  geom_boxplot()+
  ylab("Cell Count") +
  theme_bw() +
  scale_fill_manual(values =c("Total_bacteria_red" = "#FF9999", "Total bacteria_blue" =  "#56B4E9" ),
                    labels = c("Stained with DAPI", "Stained with CTC"))
 
 # percent dormancy 
 cell_stain_count_total_long |>
   ggplot(aes(Per_dormant, Species, fill = Location)) +
   geom_boxplot()
  

##################################################################
###################Microbiome ####################################
##################################################################
library(ggplot2)
library(vegan)
setwd("amphibian_dormancy/Ribbitr_data/")

#read in data
pan_meta<-read.delim("data/2022_Panama_16S_map.txt", header=T)
pan_asv<-read.delim("data/Panama_asv_table.txt", row.names = 1, header=T)

#make ASV table with negative controls only and remove from original table
pan_contam<-pan_asv[,c(which(names(pan_asv) == 'Pan27.DNA' | names(pan_asv) == 'Pan27.RNA' | names(pan_asv) == 'Pan6.DNA' | names(pan_asv) == 'Pan6.RNA' | names(pan_asv) == 'Pan4.DNA' | names(pan_asv) == 'Pan4.RNA' | names(pan_asv) == 'Pan40.DNA' | names(pan_asv) == 'Pan40.RNA'))]
pan_asv<-pan_asv[,-c(which(names(pan_asv) == 'Pan27.DNA' | names(pan_asv) == 'Pan27.RNA' | names(pan_asv) == 'Pan6.DNA' | names(pan_asv) == 'Pan6.RNA' | names(pan_asv) == 'Pan4.DNA' | names(pan_asv) == 'Pan4.RNA' | names(pan_asv) == 'Pan40.DNA' | names(pan_asv) == 'Pan40.RNA'))]

#ID which rows are contaminants
pan_contam$Sum<-rowSums(pan_contam)

#remove negative controls and filter OTU contaminated ASVs (SampleIDs = Pan27, Pan6, Pan4, Pan40)
pan_asv<-pan_asv[-which(pan_contam$Sum>0),]

#see sequencing depth
colSums(pan_asv)
length(which(colSums(pan_asv)<549))
#549 looks good, lose 11 samples (Pan38.DNA Pan41.DNA Pan44.DNA,Pan32.RNA,Pan29.RNA Pan35.RNA,Pan28.DNA,Pan24.RNA,Pan35.DNA,Pan22.DNA,Pan16.RNA)

#rarefy data
pan_rare<-rrarefy(t(pan_asv), sample=549)

#calculate PCoA based on BC similarity
pan_pcoa<-capscale(pan_rare  ~ 1, distance='bray')

#pull out x/y coordinates
pan.scores<-scores(pan_pcoa)

#grab only sample coordinates, write to data frame
pan.coords<-as.data.frame(pan.scores$sites)

#create sample names as a column
pan.coords$SampleID<-row.names(pan.coords)

#map back meta data
pan.coords<-merge(pan.coords, pan_meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(pan_pcoa$CA$eig[1]/sum(pan_pcoa$CA$eig), 3)
#30.7
100*round(pan_pcoa$CA$eig[2]/sum(pan_pcoa$CA$eig), 3)
#7

#plot PCoA
ggplot(pan.coords, aes(MDS1, MDS2, shape=Type, color=Species))+
  geom_point(aes(size=2))+
  #geom_text()+
  scale_color_manual(values = c('#f58231', '#4363d8', '#000000', "green"))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 30.7%")+
  ylab("PC2- 7%")

#create RNA/DNA only tables for calculations
library(dplyr)
pan_rna<-select(pan_asv, contains(".RNA"))
pan_dna<-select(pan_asv, contains(".DNA"))

#sort columns to be in the same order
pan_dna<-pan_dna[,order(colnames(pan_dna))]
pan_rna<-pan_rna[,order(colnames(pan_rna))]

#prune non matching samples
pan_dna<-select(pan_dna, -c(Pan22.DNA, Pan23.DNA, Pan41.DNA, Pan44.DNA, Pan25.DNA))
pan_rna<-select(pan_rna, -c(Pan24.RNA, Pan8.RNA, Pan32.RNA, Pan39.RNA, Pan8.RNA))

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
pan_dorm<-merge(pan_dorm, pan_meta, by.x='Column', by.y='SampleID')

#plot by species/site
ggplot(pan_dorm, aes(Species, Proportion, fill=region))+
  theme_bw()+
  xlab("")+
  geom_boxplot()+
  ylab("Percent active")+
  coord_flip()


