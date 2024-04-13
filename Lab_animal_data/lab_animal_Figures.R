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
library(purrr)
library(Hmisc)

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

cfu_past_fresh <- read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet =2) |>
  mutate(SampleType = "Fresh")
cfu_past_frozen <- read_excel("lab_animal_data/data/CFU_cell_counts.xlsx", sheet =6)|>
  mutate(SampleType = "Frozen")


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
 ## total_vs_viable_clean <-
  ## total_vs_viable_clean |>
  ## filter(num_spore >= 0)

total_vs_viable_clean <-
  total_vs_viable_clean |>
  mutate(fresh_vs_frozen = "Fresh")

#plotting the measure of viable versus total spore counts 
total_vs_viable_clean |>
  ggplot(aes(x= Species, y = num_spore, fill = total_vs_viable)) +
  geom_boxplot() + 
  scale_fill_manual(values =c("Total_fluor" = "#FF9999", "Viable_fluor" =  "#56B4E9" ),
                    labels = c("Viable Spore Count", "Total Spore Count"))

# plot as a percentage of viable/total
total_vs_viable |>
  ggplot(aes(x = Species, y = Per_viable, fill = Species)) +
  stat_summary(fun.data=mean_sdl, geom="bar") +
  theme_bw() +
  scale_y_continuous(labels= scales::label_percent(scale = 1, accuracy =1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  
  
# combine fresh versus frozen data for both total and viable spores
frozen_spore_total <- frozen_spore_total[-1, ]
frozen_spore_total$SampleID <- as.numeric(frozen_spore_total$SampleID)

# adds a species column based on the recorded sampleID
frozen_spore_total <- frozen_spore_total |>
  mutate(Species = case_when(
    SampleID <=5 ~ "Salamandra salamandra",
    SampleID > 5 & SampleID <= 10 ~ "Ambystoma maculatum",
    SampleID > 10 & SampleID <=15 ~ "Bombina orientalis",
    SampleID > 15 & SampleID <=20 ~ "Notophthalmus viridescens",
    SampleID > 20 & SampleID <=25 ~ "Litoria caerulea",
    SampleID > 25 & SampleID <=30 ~ "Eurycea wilderae"
  ))

# changes dataframe names to match other table for merging 
frozen_spore_total <- frozen_spore_total |>
  mutate(Fluor_reading = "Frozen")
frozen_spore_total <- setNames(frozen_spore_total, c("Well Row", "Well Col", "Content", "Frozen_fluor", 
                                                     "Type", "SampleID", "num_spore", "Species", "fresh_vs_frozen"))

# repeates steps for the viable frozen spore data
frozen_spore_viable <- frozen_spore_viable[-1, ]
frozen_spore_viable$SampleID <- as.numeric(frozen_spore_viable$SampleID)

# adds a species column based on the recorded sampleID
frozen_spore_viable <- frozen_spore_viable |>
  mutate(Species = case_when(
    SampleID <=5 ~ "Salamandra salamandra",
    SampleID > 5 & SampleID <= 10 ~ "Ambystoma maculatum",
    SampleID > 10 & SampleID <=15 ~ "Bombina orientalis",
    SampleID > 15 & SampleID <=20 ~ "Notophthalmus viridescens",
    SampleID > 20 & SampleID <=25 ~ "Litoria caerulea",
    SampleID > 25 & SampleID <=30 ~ "Eurycea wilderae"
  ))

# changes dataframe names to match other table for merging 
frozen_spore_viable <- frozen_spore_viable |>
  mutate(Type2 = "Frozen")
frozen_spore_viable <- setNames(frozen_spore_viable, c("Well Row", "Well Col", "Content", "Frozen_fluor",
                                                       "Type", "SampleID", "num_spore", "Species", "fresh_vs_frozen")) 

# create one frozen spore dataframe 
frozen_spore_combine <- bind_rows(frozen_spore_total, frozen_spore_viable)
frozen_spore_combine <- na.omit(frozen_spore_combine)

frozen_spore_combine <- frozen_spore_combine |>
  mutate(total_vs_viable = case_when(
  Type == "Swab_total" ~ "Total_fluor",
  Type == "Swab_viable" ~ "Viable_fluor",
  ))

# remove all negative spore values
frozen_spore_combine <-
  frozen_spore_combine
frozen_spore_combine$num_spore <- as.numeric(frozen_spore_combine$num_spore)
 

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

# combine fresh and frozen pasteurization dataframes
cell_count_past <- bind_rows(cfu_past_fresh, cfu_past_frozen)
View(cell_count_past)

# merge the past data with the full culture data to measure percent of total 
cell_count_past_per <- merge(Cell_culture_combined, cell_count_past) 
View(cell_count_past_per)

# graph percent of total that are spore forming bacteria 
cell_count_past_per |>
  ggplot(aes(x = Species, y = (100*(CFU_ml/`CFU/ml_2/22`)), fill = SampleType)) +
  geom_boxplot() +
  ylab("% Culturable Spore Forming Bacteria")
  



#### Standard Curve Graph ####
# create a function to apply the formula to x values
  apply_formula <- function(x, formula) {
    data.frame(x = x, y = eval(parse(text = formula), envir = list(x = x))) 
  }
  
  # Formulas to be applied
  species_formulas <- data.frame(species = c("average_curve", "Paenibacillus pabuli", "Paenibacillus tritici", 
                                             "Bacillus mycoides", "Bacillus subtilis","Paenibacillus chitinolyticus",
                                             "Bacillus_ tropicus", "Lysinibacillus fusiformis", "Virdibacillus arri"),
  formulas <- c("0.272*x+22798", "0.277*x+24384", "0.263*x+23166", "0.254*x+21567", 
                "0.277*x+17471", "0.291*x+25604", "0.291*x+27622", "0.291*x+21181", 
                "0.272 * x + 23898")
  )
  
 # Generate random x values
  x_values <- runif(100, 0, 15000)
  
  # Apply formulas and create a data frame
  result_df <- map2_df(species_formulas$formula, species_formulas$species, ~ apply_formula(x_values, .x) |>
                         mutate(species = .y, formula = .x))

 #vizualize result
  View(result_df)
  
 # make a plot with error bars using the newly created data
  result_df |>
    ggplot(aes(x=x, y=y, color = species, group = species)) +
    geom_line()+
    geom_smooth(method="lm") +
    theme_minimal()
  










