#----------------------------
#
# Andrew Miller-Klugman
# Ribbitr Data Figures
#
#----------------------------

library(dyply)
library(ggplot2)
library(readxl)
library(tidyr)

# load in the data
cell_counts <- read_xlsx("Ribbitr_data/data/RIBBiTR_CFUs_cell_counts.xlsx")
TbCl_data <- read_xlsx("Ribbitr_data/data/RIBBiTR_TbCl_data.xlsx")

#### Cell Culture Count Data ####
# boxplot of bacterial CFUs by species per location
cell_counts |>
  ggplot(aes(Species...5, CFU_bact...25, fill = Location...2)) +
  geom_boxplot()

#boxplot of pasteurized CFUs
cell_counts |>
  ggplot( aes(Species...5, CFU_past, fill = Location...2)) +
  geom_boxplot()
  
# Boxplot of %spore forming bacteria 
cell_counts |>
  ggplot(aes(Species...5, Percent_spore, fill = Location...2)) +
  geom_boxplot(outlier.shape = NA, outlier.colour = NA) +
  scale_y_continuous(limits = c(0,10))

# Percent_spore by month
cell_counts |>
  ggplot(aes(Species...5, Percent_spore, fill = Month_sampled...6)) +
  geom_boxplot(outlier.shape = NA, outlier.colour = NA) +
  scale_y_continuous(limits = c(0,10))

# s16 copies by species
cell_counts |>
  ggplot(aes(Species...5, S16_copies, fill = Location...2)) +
  geom_boxplot()


#### TbCl Data

#percent viable spore
TbCl_data |>
  ggplot(aes(Species, Per_viable_spore, fill = Location)) +
  geom_boxplot()
  
