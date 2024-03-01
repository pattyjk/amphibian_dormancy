#----------------------------
#
# Andrew Miller-Klugman
# Lab Animal Data Figures
#
#----------------------------

# Load in Libraries
library(readr)
library(readxl)
library(ggplot2)
library(tidyr)
library(data.table)
library(dplyr)
library(agricolae)

# Load in the needed data
#TbCl Data 
viable_spores <-
  read_excel("lab_animal_data/data/TbCl_swabs_viable_spores.xlsx")

total_spore <-
  read_excel("lab_animal_data/data/tbcl_total_spores_swab.xlsx")

frozen_spore_viable <-
  read_excel("lab_animal_data/data/TbCl_assay_post_freeze.xlsx", sheet = 4)

frozen_spore_total <-
  read_excel("lab_animal_data/data/TbCl_assay_post_freeze.xlsx", sheet = 2)

total_vs_viable <- read.delim("lab_animal_data/data/per_viable_spores.txt", sep = '\t', header=T)

#cell stain data 
cell_stains <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet = 7)

#cell culture data 
cell_culture_fresh <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet = 3)

cell_culture_frozen <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet =5)



#### Cell Stain Figures ####

#pivot data to longer to place on a singular plot 
cell_stain_long <-
  pivot_longer(cell_stains, cols = c('Total bacteria_blue', 'Total_bacteria_red'), names_to = "stain", values_to = "value")

# create a barplot that displays both stain types on a single axis based on species
ctc_dapi <-
  ggplot(cell_stain_long, aes(x = value, y = Species, fill = stain)) +
           geom_boxplot()+
  ylab("Cell Count") +
  theme_bw() +
  scale_fill_manual(values =c("Total_bacteria_red" = "#FF9999", "Total bacteria_blue" =  "#56B4E9" ),
                    labels = c("Stained with DAPI", "Stained with CTC"))

# Stats 
set.seed(12345)
# run anova to compare effects of stain type and species 
anova_result_stains <- aov(value ~ stain * Species, data = cell_stain_long)
summary(anova_result_stains)

# TukeyHSD to view differences by species 
tukey_result_stain <- TukeyHSD(anova_result_stains, "Species", group = TRUE)
print(tukey_result_stain)


#### TbCL Figures ####

# Remove mucosome data
total_vs_viable <-
  total_vs_viable |>
  filter(Type1 != 'Mucosome')

# pivot datato longer to place on singular plot
total_vs_viable_long <-
  pivot_longer(total_vs_viable, cols = c('Total_fluor', 'Viable_fluor'), names_to = "total_vs_viable", values_to = "Fluor_reading")

#adds a new column that determines the spore count from the fluorescnese reading data based on the average standard curve
total_vs_viable_clean <-
  total_vs_viable_long |>
  mutate(num_spore = (Fluor_reading-22798)/0.272)

#remove the negative values from the spore counts
total_vs_viable_clean <-
  total_vs_viable_clean |>
  filter(num_spore >= 0)


#plotting the measure of viable versus total spore counts 
total_vs_viable_clean |>
  ggplot(aes(x= Species, y = num_spore, fill = total_vs_viable)) +
  geom_boxplot() + 
  scale_fill_manual(values =c("Total_fluor" = "#FF9999", "Viable_fluor" =  "#56B4E9" ),
                    labels = c("Viable Spore Count", "Total Spore Count"))
  

# Stats 
anova_result_tbcl <- aov(num_spore ~ total_vs_viable * Species, data = total_vs_viable_clean)
summary(anova_result_tbcl)

# TukeyHSD to view differences by species 
tukey_result_tbcl <- TukeyHSD(anova_result_tbcl, "Species", group = TRUE)
print(tukey_result_tbcl)


#### Cell Culture Figures #####

# adds a new column denoting the samples as fresh or frozen for combination into a single table

cell_culture_frozen<- cell_culture_frozen |>
  mutate(SampleType = "Frozen")

cell_culture_fresh<- cell_culture_fresh |>
  mutate(SampleType = "Fresh")

# Cleaning Data table
cell_culture_fresh <- cell_culture_fresh |>
  select(`CFU/ml_2/22`, `CFUs_2/22`, Sample_ID, SampleType, Species, Common_name)

cell_culture_frozen <- cell_culture_frozen |>
  select(`CFU/ml_2/22`,  Colonies, Sample_ID, SampleType, Species, Common_name) |>
  rename('CFUs_2/22' = Colonies)

#Combine the two data tables into one for making figures
Cell_culture_combined <- 
  bind_rows(cell_culture_fresh, cell_culture_frozen)

#Figure showing CFUS/ml for both fresh and frozen samples per species

Cell_culture_combined |>
  ggplot( aes(`CFU/ml_2/22`, Species , fill = SampleType)) +
  geom_boxplot() +
  theme_bw()


#### Standard Curve Graph ####

#create functions for each standard curve
Virdibacillus_arri <- function (x) {
   return(0.272 * x + 23898)
}

Lysinibacillus_fusiformis <- function(x) {
  return(0.291*x+21181)
}

Bacillus_tropicus <- function(x) {
  return(0.291*x+27622)
}

Paenibacillus_chitinolyticus <- function(x) {
  return(0.291*x+25604)
}

Bacillus_subtilis <- function(x) {
  return(0.277*x+17471)
}

Bacillus_mycoides <- function(x) {
  return(0.254*x+21567)
}

Paenibacillus_tritici <- function(x) {
  return(0.263*x+23166)
}

Paenibacillus_pabuli <- function(x) {
  return(0.277*x+24384)
}

average_curve <- function(x) {
  return(0.272*x+22798)
}

x_values <- seq(0, 15000, length.out = 100)
data <- data.frame(x = x_values)

#plot the standard curves on the graph 
  ggplot(data, aes(x=x)) +
  geom_line(aes(y = Virdibacillus_arri(x), color = "Virdibacillus arri"), linewidth = 1) +
    geom_line(aes(y =Lysinibacillus_fusiformis(x), color = "ysinibacillus fusiformis"), linewidth = 1) +
    geom_line(aes(y = Bacillus_tropicus(x), color = "Bacillus tropicus"), linewidth = 1) +
    geom_line(aes(y = Paenibacillus_chitinolyticus(x), color = "Paenibacillus chitinolyticus"), linewidth = 1) +
    geom_line(aes(y = Bacillus_subtilis(x), color = "Bacillus subtilis"), linewidth = 1) +
    geom_line(aes(y = Bacillus_mycoides(x), color = "Bacillus mycoides"), linewidth = 1) +
    geom_line(aes(y = Paenibacillus_tritici(x), color = "Paenibacillus tritici"), linewidth = 1) +
    geom_line(aes(y = Paenibacillus_pabuli(x), color = "Paenibacillus pabuli"), linewidth = 1) + 
    geom_line(aes(y = average_curve(x), color = "Average curve"), linewidth = 1) +
  labs(title = "Standard Curves", x = "Spore Count", y = "Fluorescence") +
  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(family = 'serif', face = 'bold', size = 12))
  
  














