library(geodata)
library(terra)

# 2. Load your coordinate file
# Assumes 'ribbitr_latlong.txt' is in your current working directory
coords_df <- read.table("~/Documents/GitHub/amphibian_dormancy/ribbitr_latlong.txt", header = TRUE, sep = "\t")

#Download WorldClim bioclimatic data, takes time to grab
# res = 0.5 (highest resolution, ~1km), 2.5, 5, or 10 minutes of a degree
# var = 'bio' downloads all 19 bioclimatic variables
bioclim_data <- worldclim_global(var = "bio", res = 0.5, path = "data/")

#Extract climate data for your coordinates
extracted_values <- extract(bioclim_data, coords_df[, c("site_longitude", "site_latitude")])

#Combine original coordinates with the new climate data
final_data <- cbind(coords_df, extracted_values)

#add meaningful column names
# Descriptive names for the 19 Bioclimatic variables
bio_colnames <- c(
  "lat", 'long', 'id',
  "annual_mean_temp",          # BIO1
  "mean_diurnal_range",        # BIO2
  "isothermality",             # BIO3
  "temp_seasonality",          # BIO4
  "max_temp_warmest_month",    # BIO5
  "min_temp_coldest_month",    # BIO6
  "temp_annual_range",         # BIO7
  "mean_temp_wettest_quarter", # BIO8
  "mean_temp_driest_quarter",  # BIO9
  "mean_temp_warmest_quarter", # BIO10
  "mean_temp_coldest_quarter", # BIO11
  "annual_precip",             # BIO12
  "precip_wettest_month",      # BIO13
  "precip_driest_month",       # BIO14
  "precip_seasonality",        # BIO15
  "precip_wettest_quarter",    # BIO16
  "precip_driest_quarter",     # BIO17
  "precip_warmest_quarter",    # BIO18
  "precip_coldest_quarter"     # BIO19
)

names(final_data)<-bio_colnames

#write to file for safe keeping
write.table(final_data, '~/Documents/GitHub/amphibian_dormancy/ribbitr_climate_data.txt', sep='\t', quote=F, row.names=F)
