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

#### Cell Stain Figures ####

ctc_dapi <-
  ggplot(cell_stains, mapping = aes(x = Species, y = cell_stains$`Total bacteria_blue`) |>
           geom_boxplot()


  

