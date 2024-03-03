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
  
