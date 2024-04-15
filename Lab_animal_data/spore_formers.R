#'------------------------------------------------------------
#'             Andrew Miller-Klugman
#'        Spore forming Genera Abundance
#'              Started: 2024-03-23
#'-----------------------------------------------------------

# read in the file of the spore forming genera 
spore_forming_genera <- readLines("spore_forming_genera.txt")
# load in the asv files
source("Lab_animal_data/Microbiome_data.R")

taxonomy <- read_tsv("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Lab_animal_data/data/taxonomy.tsv")

# add a total row to the taxon df that contains the total for each column
sum_row <- colSums(asv_table)

# Convert the sums to a dataframe
sum_row <- as.data.frame(t(sum_row))
rownames(sum_row) <- "Total Bacteria Count"

asv_spore <- rbind(asv_table, sum_row)

# cleaning the asvtable data to allow for merging
asv_spore <- rownames_to_column(asv_table, var = "Feature ID")
taxon <- merge(asv_spore, taxonomy, by = "Feature ID")

# filters for spore forming genera based on whats listed in the textfile
spore_formers1 <- taxon[grep(paste(spore_forming_genera, collapse="|"), taxon$Taxon), ]
spore_formers <- select(spore_formers1, -c(Taxon, Confidence))
rownames(spore_formers) <- c("053480b13b8379901eb13a509889164a", "14d168c89e381de6fc3091837a557e48", "d0958c80b3756a5cf6cf4e71c5693e64")
#remove feature ID column for combining with total data 
spore_formers_fixed <- select(spore_formers, -c(`Feature ID`))
#add the total sums of all bacteria counts to the spore formers table
spore_formers_final <- rbind(spore_formers_fixed, sum_row)
# calculate percentages of spore formers
per_spore <- colSums(spore_formers_final[1:3,]) / spore_formers_final[4,] *100
rownames(per_spore) <- "Percent Spore Forming Bacteria"
# combine the percent spore formers with the main table
per_spore_form <- rbind(spore_formers_final, per_spore)

# create a plot showing percentage of spore forming bacteria 
# add metadata to the spore table
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

# change all number values into type numeric, once again because R is dumb 
numeric_per_spore <- mutate_at(per_spore_filter, .vars= c("053480b13b8379901eb13a509889164a", 
                                                         "14d168c89e381de6fc3091837a557e48", "d0958c80b3756a5cf6cf4e71c5693e64", 
                                                         "Percent Spore Forming Bacteria", "Total Bacteria Count"),
                               as.numeric)

# plot the percentage of spore forming bacteria per species
dormant_16s <- numeric_per_spore |>
ggplot(aes(x = Species, y = `Percent Spore Forming Bacteria`, fill = Species)) +
  geom_boxplot()

# output a table of spore forming bacteria that are present in the species
spore_forming_bacteria <- select(spore_formers1, c("Feature ID", "Taxon", "Confidence"))
write.table(spore_forming_bacteria, "Lab_animal_data/data/Lab_animal_16s_spore_formers.txt", sep='\t', quote=F, row.names=F)
write.table(numeric_per_spore, "Lab_animal_data/data/sample_data_meta.csv", sep = ",", quote=F, row.names=F)

dat <- asv_spore |>
  pivot_longer(-`Feature ID`, names_to = "SampleID", values_to = "count")

# merge taxonomy
dat <- dat |> 
  left_join(merged_taxonomy, by = "Feature ID")

sample_meta <- read_csv("Lab_animal_data/data/sample_data_meta.csv")
meta_dat <- dat |>
  left_join(sample_meta, by = "SampleID")

meta_dat <- meta_dat[!(meta_dat$Species == "" | 
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
    Species == "Bombina orientalis" ~ "B. orientalis",
    Species == "Notophthalmus viridescens" ~ "N. viridescens",
    Species == "Salamandra salamandra" ~ "S. salamandra",
    Species == "Ambystoma maculatum" ~ "A. maculatum",
    Species == "Eurycea wilderae" ~ "E. wilderae",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))


 sample_meta_clean |>
  ggplot(aes(x = `Species Name`, y = count)) +
  geom_bar(aes(fill = genus), stat = "identity", position = "fill") +
  scale_y_continuous(name = "Relative Abundance") +
  labs(fill = "Bacteria Family") +
   theme_minimal() +
  theme(text = element_text(size = 32, family = "Arial")) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(face = "italic")) +
   scale_fill_manual(values = c("#B79F00", "#00BA38", "#619CFF"))
 
 # plot entire bacterial composition 
 
 # Function to extract family name
 extract_phylum <- function(taxonomy_chain) {
   # Split the taxonomy chain into components
   components <- strsplit(taxonomy_chain, "; ")[[1]]
   # Find the component containing the family name
   family_component <- grep("^p__", components)
   # Extract the family name
   if (length(family_component) > 0) {
     family <- gsub("^p__", "", components[family_component])
     family <- gsub("_.*", "", family)
     return(family)
   } else {
     return(NA)  # If family name is not found, return NA
   }
 }
 
 meta_dat$genus <- sapply(meta_dat$Taxon, extract_phylum)
 
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
     Species == "Bombina orientalis" ~ "B. orientalis",
     Species == "Notophthalmus viridescens" ~ "N. viridescens",
     Species == "Salamandra salamandra" ~ "S. salamandra",
     Species == "Ambystoma maculatum" ~ "A. maculatum",
     Species == "Eurycea wilderae" ~ "E. wilderae",
     TRUE ~ NA_character_  # Handle other cases if needed
   ))
 
 meta_dat_1 |>
   ggplot(aes(x = `Species Name`, y = count)) +
   geom_bar(aes(fill = `genus`), stat = "identity", position = "fill") +
   scale_y_continuous(name = "Relative Abundance") +
   labs(fill = "Bacteria Family") +
   theme_minimal() +
   theme(text = element_text(size = 32, family = "serif")) +
   theme(strip.background = element_blank(),
         axis.text.x = element_text(face = "italic"))
