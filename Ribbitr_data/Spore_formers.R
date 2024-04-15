#` ------------------------------------------------
#'           Andrew Miller-Klugman
#'         RIBBiTR 16s Spore Formers
#'            Started: 2024-03-30
#' ------------------------------------------------

source("Ribbitr_data/RIBBiTR_microbiome.R")

taxonomy1 <- read_tsv("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Ribbitr_data/data/exported-feature-table/taxonomy1/taxonomy.tsv")
taxonomy2 <- read_tsv("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Ribbitr_data/data/exported-feature-table/taxonomy2/taxonomy.tsv")
merged_taxonomy <- rbind(taxonomy1, taxonomy2, by = "Feature ID")

# reads in txt file of spore forming genera
spore_forming_genera <- readLines("spore_forming_genera.txt")

# add a total row to the taxon df that contains the total for each column
sum_row <- colSums(asv_table)

# Convert the sums to a dataframe
sum_row <- as.data.frame(t(sum_row))
rownames(sum_row) <- "Total Bacteria Count"

asv_spore <- rbind(asv_table, sum_row)
sum_row_f <- rownames_to_column(sum_row, var = "Feature ID")

# cleaning the asvtable data to allow for merging
asv_spore <- rownames_to_column(asv_spore, var = "Feature ID")
taxon <- left_join(asv_spore, merged_taxonomy, by = "Feature ID")

# filters for spore forming genera based on whats listed in the textfile
spore_formers1 <- taxon[grep(paste(spore_forming_genera, collapse="|"), taxon$Taxon), ]

# adding total counts to the dataframe 
spore_formers <- select(spore_formers1, -c(Taxon, Confidence))
# add the total sums of all bacteria counts to the spore formers table
spore_formers_final <- rbind(spore_formers, sum_row_f)
spore_formers_final <- remove_rownames(spore_formers_final)
spore_formers_final <- distinct(spore_formers_final)

# add percentage of spore forming bacteria basd on total counts
spore_percentage <- column_to_rownames(spore_formers_final, var = "Feature ID")
per_spore <- colSums(spore_percentage[1:95,]) / spore_percentage[96,] *100
rownames(per_spore) <- "Percent Spore Forming Bacteria"

# combine the percent spore formers with the main table
per_spore_form <- rbind(spore_percentage, per_spore)

# add metadata to the spore_forming percentage data
# Transpose the dataset
transposed_data <- t(meta)

# Convert the transposed data to a dataframe
transposed_data <- as.data.frame(transposed_data)

# Set the column names to be the first row values
colnames(transposed_data) <- transposed_data[1, ]

# Remove the first row
transposed_data <- transposed_data[-1, ]

# Print the transposed dataset
print(transposed_data)

common_columns <- intersect(names(per_spore_form), names(transposed_data))

# Convert row names into a column
transposed_data$Row.names <- row.names(transposed_data)
per_spore_form$Row.names <- row.names(per_spore_form)
per_spore_meta <- merge(per_spore_form, transposed_data, by = c("Row.names", common_columns), all = TRUE)
rownames(per_spore_meta) <- per_spore_meta$Row.names
per_spore_meta$Row.names <- NULL  # Remove the auxiliary column

# transpose format for graphing
per_spore_transpose <- t(per_spore_meta)
# Convert the transposed data to a dataframe
per_spore_transpose <- as.data.frame(per_spore_transpose)
# make the rownames into a column because R is dumb sometimes
per_spore_transpose <- rownames_to_column(per_spore_transpose, var = "SampleID")
# remove NA values
per_spore_filter <- na.omit(per_spore_transpose)

# change all numbers to a numeric case and remove all empty values
numeric_per_spore <- mutate_if(per_spore_filter, ~all(grepl("^\\d+\\.?\\d*$", .)), as.numeric)
numeric_per_spore_clean <- numeric_per_spore[!(numeric_per_spore$Location == "" |numeric_per_spore$Species == "" | 
                                                 numeric_per_spore$Species == "N/A" | numeric_per_spore$Species =="Negative Control" |
                                                 numeric_per_spore$Species =="Negative control"), ]

# plot the percentage of spore forming bacteria per species
per_dormant_16s <- numeric_per_spore_clean |>
  ggplot(aes(x = Species, y = `Percent Spore Forming Bacteria`, fill = Species)) +
  geom_boxplot()

numeric_per_spore_clean |>
  ggplot(aes(x = Location, y = `Percent Spore Forming Bacteria`, fill = Species)) +
  geom_boxplot()

# merge the tacon data back with the percent spore formers
final_metadata <- rownames_to_column(per_spore_filter, var = "Feature ID")

write.table(numeric_per_spore_clean, "Ribbitr_data/data/RIBBiTR_16s_spore_formers.csv", sep=',', quote=F, row.names=F)
write.table(sample_meta_clean, "Ribbitr_data/data/sample_meta_clean.csv", sep =',', quote=F, row.names=F)
write.table(meta_dat, "Ribbitr_data/data/sample_meta.csv", sep =',', quote=F, row.names=F)

# relative abundance of spore forming bacteria per species 
dat <- asv_spore |>
  pivot_longer(-`Feature ID`, names_to = "SampleID", values_to = "count")

# merge taxonomy
dat <- dat |> 
left_join(merged_taxonomy, by = "Feature ID")

sample_meta <- read_csv("Ribbitr_data/data/sample_meta_data.csv")
meta_dat <- dat |>
  left_join(sample_meta, by = "SampleID")

# remove all NA values 
meta_dat <- meta_dat[!(meta_dat$Location == "" |meta_dat$Species == "" | 
                         meta_dat$Species == "N/A" | meta_dat$Species =="Negative Control" |
                         meta_dat$Species =="Negative control" | is.na(meta_dat$Species)), ]

sample_meta_clean <- subset(meta_dat, grepl(paste(spore_forming_genera, collapse="|"), meta_dat$Taxon))

# Function to extract family name
extract_family <- function(taxonomy_chain) {
  # Split the taxonomy chain into components
  components <- strsplit(taxonomy_chain, "; ")[[1]]
  # Find the component containing the family name
  family_component <- grep("^f__", components)
  # Extract the family name
  if (length(family_component) > 0) {
    family <- gsub("^f__", "", components[family_component])
    family <- gsub("_.*", "", family)
    return(family)
  } else {
    return(NA)  # If family name is not found, return NA
  }
}

# Apply the function to each row in the dataframe
sample_meta_clean$genus <- sapply(sample_meta_clean$Taxon, extract_family)

# add scientific name for plotting ease
sample_meta_clean <- sample_meta_clean |>
  mutate(`Species Name` = case_when(
    Species == "Rana pipiens" ~ "R. pipiens",
    Species == "Rana catesbeiana" ~ "R. catesneiana",
    Species == "Hylodes_phyllodes" ~ "H. phyllodes",
    Species == "Ischnocnema_henselii" ~ "I. Henselii",
    Species == "Colostethus panamansis" ~ "C. panamansis",
    Species == "Lithobates warzsewitschii" ~ "L. warszewitschii",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))


# plotting the relative abundance of each spore forming bacteria 
sample_meta_clean |>
  ggplot(aes(x = `Species Name`, y = count)) +
  geom_bar(aes(fill = genus), stat = "identity", position = "fill") +
  facet_grid(~ Location, scales = "free_x", space = "free_x") +
  scale_y_continuous(name = "Relative Abundance") +
  labs(fill = "Bacteria Family") +
  theme(text = element_text(size = 14, family = "serif", face = "Bold")) +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(face = "italic"))
  

# Function to extract family name
extract_order <- function(taxonomy_chain) {
  # Split the taxonomy chain into components
  components <- strsplit(taxonomy_chain, "; ")[[1]]
  # Find the component containing the family name
  family_component <- grep("^o__", components)
  # Extract the family name
  if (length(family_component) > 0) {
    family <- gsub("^o__", "", components[family_component])
    family <- gsub("_.*", "", family)
    return(family)
  } else {
    return(NA)  # If family name is not found, return NA
  }
}

meta_dat$genus <- sapply(meta_dat$Taxon, extract_order)

meta_dat_1 <- meta_dat |>
  mutate(`spore_form` = case_when(
    str_detect(Taxon, "Actinmycetaceae") ~ "Spore Former",
    str_detect(Taxon, "Bacillaceae") ~ "Spore Former",
    str_detect(Taxon, "Clostridiacease") ~ "Spore Former",
    str_detect(Taxon, "Desulfitobacteriaceae") ~ "Spore Former",
    str_detect(Taxon, "Paenibacillaceae") ~ "Spore Former",
    str_detect(Taxon, "Syntrophomonadaceae") ~ "Spore Former",
    TRUE ~ "Not Spore Former"  # Handle other cases if needed
  )
  )

 meta_dat$genus <- sapply(meta_dat$Taxon, extract_order)
 
 meta_dat_1 <- meta_dat |>
   mutate(`spore_form` = case_when(
     str_detect(Taxon, "Actinmycetaceae") ~ "Spore Former",
     str_detect(Taxon, "Bacillaceae") ~ "Spore Former",
     str_detect(Taxon, "Clostridiacease") ~ "Spore Former",
     str_detect(Taxon, "Desulfitobacteriaceae") ~ "Spore Former",
     str_detect(Taxon, "Paenibacillaceae") ~ "Spore Former",
     str_detect(Taxon, "Syntrophomonadaceae") ~ "Spore Former",
     TRUE ~ "Not Spore Former"  # Handle other cases if needed
   )
   )
 
 # add scientific name for plotting ease
 meta_dat_1 <- meta_dat_1 |>
   mutate(`Species Name` = case_when(
     Species == "Rana pipiens" ~ "R. pipiens",
    Species == "Rana catesbeiana" ~ "R. catesneiana",
    Species == "Hylodes_phyllodes" ~ "H. phyllodes",
    Species == "Ischnocnema_henselii" ~ "I. Henselii",
    Species == "Colostethus panamansis" ~ "C. panamansis",
    Species == "Lithobates warzsewitschii" ~ "L. warszewitschii",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))


 meta_dat_1 |>
   ggplot(aes(x = `Species Name`, y = count)) +
   geom_bar(aes(fill = `spore_form`), stat = "identity", position = "fill") +
   scale_y_continuous(name = "Relative Abundance") +
   labs(fill = "Bacteria Family") +
   theme_minimal() +
   theme(text = element_text(size = 32, family = "Arial")) +
   theme(strip.background = element_blank(),
         axis.text.x = element_text(face = "italic"))

#### FINAL BARPLOT ####

# entire bacterial community 
final_ribbitr <- read.csv("Ribbitr_data/data/sample_meta_clean.csv")
total_abundance <- meta_dat |>
  mutate(total_abundance = sum(count))

final_ribbitr |>
  ggplot(aes(x = Species.Name, y = Relative.Abundance*100)) +
  geom_bar(aes(fill = `genus`), stat = "identity", position = "fill") +
  scale_y_continuous(name = "Relative Abundance") +
  labs(fill = "Bacteria Family") +
  theme_minimal() +
  theme(text = element_text(size = 28, family = "serif")) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(face = "italic"))


meta_dat_1 |>
  ggplot(aes(x = `Species Name`, y = count)) +
  geom_bar(aes(fill = `spore_form`), stat = "identity", position = "fill") +
  scale_y_continuous(name = "Relative Abundance") +
  labs(fill = "Bacteria Family") +
  theme_minimal() +
  theme(text = element_text(size = 32, family = "Arial")) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(face = "italic"))

