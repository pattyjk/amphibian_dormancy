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


#### Example ####

#calculate correlations

#read in data
data<-read.delim("amphibian_dormancy/sample_data.txt", header=T)
str(dat)

##write a function to do the calculation (SPearman correlation)
#parts to change: 
#data- whatever your data frame is named, 
#the 'fresh_CFU' and 'Frozen_counts' to whatever column you're naming
#under part 2, change the species to whatever you're comparing

calculate_correlation <- function(factor_level, data) {
  subset_data <- subset(data, Species == factor_level)
  correlation_test <- cor.test(subset_data$Frozen_counts, subset_data$Fresh_CFU, method='spearman')
  return(correlation_test)
}

#Calculate correlation with p-value for each level of the factor
factor_levels <- unique(data$Species)
correlation_results <- lapply(factor_levels, function(level) {
  calculate_correlation(level, data)
})

#view results
correlation_results



