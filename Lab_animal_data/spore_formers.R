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
numeric_per_spore |>
ggplot(aes(x = Species, y = `Percent Spore Forming Bacteria`, fill = Species)) +
  geom_boxplot()

# output a table of spore forming bacteria that are present in the species
spore_forming_bacteria <- select(spore_formers1, c("Feature ID", "Taxon", "Confidence"))
write.table(spore_forming_bacteria, "Lab_animal_data/data/Lab_animal_16s_spore_formers.txt", sep='\t', quote=F, row.names=F)

