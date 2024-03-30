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

# cleaning the asvtable data to allow for merging
asv_spore <- rownames_to_column(asv_table, var = "Feature ID")
taxon <- merge(asv_spore, taxonomy, by = "Feature ID")

# filters for spore forming genera based on whats listed in the textfile
spore_formers <- taxon[grep(paste(spore_forming_genera, collapse="|"), taxon$Taxon), ]

# add a total row to the taxon df that contains the total for each column
total_counts <- colSums(taxon[, -1])
