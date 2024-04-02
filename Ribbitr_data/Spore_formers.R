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

numeric_per_spore <- mutate_if(per_spore_filter, ~all(grepl("^\\d+\\.?\\d*$", .)), as.numeric)

# plot the percentage of spore forming bacteria per species
numeric_per_spore |>
  ggplot(aes(x = Species, y = `Percent Spore Forming Bacteria`, fill = Species)) +
  geom_boxplot()

# merge the tacon data back with the percent spore formers

final_metadata <- rownames_to_column(per_spore_filter, var = "Feature ID")
RIBBiTR_Microbiome <- left_join(final_metadata, merged_taxonomy, by = "Feature ID")


#### TEST ####
test_df <- t(spore_formers1)
test_df <- as.data.frame((test_df))
test_df <- rownames_to_column(test_df, var = "Sample ID")
colnames(test_df) <- test_df[1, ]

# Remove the first row
test_df <- test_df[-1, ]


rows_to_transpose <- 464:465  # Example: transposing the first two rows

# Transpose selected rows
transposed_rows <- t(test_df[rows_to_transpose, ])

# Convert transposed rows to dataframe
transposed_test_df <- as.data.frame(transposed_rows)
transposed_Test_df <- rownames_to_column(transposed_test_df, var = )

# Remove transposed rows from original dataframe
test_df <- test_df[-rows_to_transpose, ]

# Combine transposed dataframe with original dataframe
result_df <- rbind(test_df, transposed_test_df)

# Print the result
print(result_df)