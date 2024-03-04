# --------------------------------------------------
#
#   Andrew Miller-Klugman
#   Correlational Table
#   Started: 2023-02-01
#
# --------------------------------------------------

# Load in Libraries
library(readr)
library(readxl)
library(ggplot2)
library(tidyr)
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

#cell stain data 
cell_stains <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet = 7)

#cell culture data 
cell_culture_fresh <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet = 3)

cell_culture_frozen <-
  read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet =5)



#### format data correctly ####

cell_culture_frozen<- cell_culture_frozen |>
  rename(Frozen_counts = 'CFU/ml_2/22')

cell_culture_fresh<- cell_culture_fresh |>
  rename(Fresh_counts = 'CFU/ml_2/22')

# Cleaning Data table
cell_culture_fresh <- cell_culture_fresh |>
  select(Fresh_counts, `CFUs_2/22`, Sample_ID, Species, Common_name)

cell_culture_frozen <- cell_culture_frozen |>
  select(Frozen_counts, Sample_ID, Species, Common_name)

#Combine the two data tables into one for making figures
Cell_culture_combined <- 
  merge(cell_culture_fresh, cell_culture_frozen, by = c("Species", "Sample_ID"))

Cell_culture_combined <-
  Cell_culture_combined |>
  select(Fresh_counts, Frozen_counts, Sample_ID, Species)
#### Example ####

#calculate correlations

##write a function to do the calculation (SPearman correlation)
#parts to change: 
#data- whatever your data frame is named, 
#the 'fresh_CFU' and 'Frozen_counts' to whatever column you're naming
#under part 2, change the species to whatever you're comparing


#### CELL CULTURE DATA ####
calculate_correlation <- function(factor_level, Cell_culture_combined) {
  subset_data <- subset(Cell_culture_combined, Species == factor_level)
  correlation_test <- cor.test(subset_data$Frozen_counts, subset_data$Fresh_counts, method='spearman')
  return(correlation_test)
}

#Calculate correlation with p-value for each level of the factor
factor_levels <- unique(Cell_culture_combined$Species)
correlation_results <- lapply(factor_levels, function(level) {
  calculate_correlation(level, Cell_culture_combined)
})

#view results
correlation_results

#### TBCL ASSAY DATA ####
#Calculate correlation with p-value for each level of the factor
factor_levels <- unique(frozen_fresh_tbcl$Species)
correlation_results <- lapply(factor_levels, function(level) {
  calculate_correlation(level, frozen_fresh_tbcl)
})



