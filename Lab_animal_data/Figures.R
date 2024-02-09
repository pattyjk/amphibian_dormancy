#----------------------------
#
# Andrew Miller-Klugman
# Fresh versus Frozen Data Figures
#
#----------------------------

# Load in Libraries
library(readr)
library(readxl)
library(ggplot2)
library(tidyr)
library(data.table)
library(dplyr)

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
  geom_boxplot()


# Remove the Mucosome data
tbcl_spores_filtered <- subset(total_vs_viable, !grepl("Mucosome", Type1, ignore.case = TRUE))


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





