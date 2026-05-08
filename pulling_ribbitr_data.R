#Pulling data (mucosome, environmental, amp) from RIBBiTR database for microbiome samples
library(ribbitrrr)
library(tidyverse)
library(dbplyr)
library(RPostgres) 
library(DBI)

#only if need to reconfigure connection
usethis::edit_r_environ()

#connect to database
dbcon = hopToDB("ribbitr")

# load table "all_tables" from schema "public"
mdt <- tbl(dbcon, Id("public", "all_tables")) %>%
  collect()

# load table "all_columns" from schema "public", filtering to schema "survey_data"
mdc <- tbl(dbcon, Id("public", "all_columns")) %>%
  filter(table_schema == "survey_data") %>%
  collect()

db_capture = tbl(dbcon, Id("survey_data", "capture"))
db_survey = tbl(dbcon, Id("survey_data", "survey"))
db_visit = tbl(dbcon, Id("survey_data", "visit"))
db_site = tbl(dbcon, Id("survey_data", "site"))
db_region = tbl(dbcon, Id("survey_data", "region"))
db_country = tbl(dbcon, Id("survey_data", "country"))
db_survey = tbl(dbcon, Id("survey_data", "survey"))
db_visit = tbl(dbcon, Id("survey_data", "visit"))
db_site = tbl(dbcon, Id("survey_data", "site"))
db_region = tbl(dbcon, Id("survey_data", "region"))
db_country = tbl(dbcon, Id("survey_data", "country"))
db_envir = tbl(dbcon, Id("survey_data", "environmental"))

# sample tables
db_sample = tbl(dbcon, Id("survey_data", "sample"))
db_bdqpcr = tbl(dbcon, Id("survey_data", "bd_qpcr_results"))
db_amp_total = tbl(dbcon, Id("survey_data", "amp_total"))
db_amp_gia = tbl(dbcon, Id("survey_data", "amp_gia"))
db_amp_peak = tbl(dbcon, Id("survey_data", "amp_maldi_peak"))
db_amp_intensity = tbl(dbcon, Id("survey_data", "amp_maldi_intensity"))
db_mucosome_gia = tbl(dbcon, Id("survey_data", "mucosome_gia"))

# Build data query
# all samples
capture_all = db_capture %>%
  left_join(db_survey, by = "survey_id") %>%
  left_join(db_visit, by = "visit_id") %>%
  left_join(db_site, by = "site_id") %>%
  left_join(db_region, by = "region_id") %>%
  left_join(db_country, by = "country_id") %>%
  left_join(db_envir, by = "survey_id")

# pull out specific sample types, joining with desired results tables (where present)
sample_mb = db_sample %>%
  filter(sample_type == "microbiome")

result_bd = db_sample %>%
  filter(sample_type == "bd") %>%
  inner_join(db_bdqpcr %>%
               filter(!is.na(bd_detected)) %>%
               group_by(sample_id, sample_name_bd) %>%  # for some sets of results we we have replicates, and need to handle these appropriately for your needs
               summarise(bd_replicates = n(),
                         bd_detected = any(bd_detected, na.rm = TRUE),
                         bd_mean_its1_copies_per_swab = mean(bd_its1_copies_per_swab, na.rm = TRUE),
                         .groups = "drop"), by = "sample_id")


# join with whichever amp results tables are of interest (total peptides, GIA, maldi peak, maldi intensity). Note that GIA and MALDI data need aggregation
result_amp = db_sample %>%
  filter(sample_type == "amp") %>%
  inner_join(db_amp_total, by = "sample_id")

# join all samples together by capture id
sample_cols_drop = c("sample_name",
                     "sample_type",
                     "negative_control",
                     "negative_control_group_id",
                     "sample_name_conflict")

capture_sample_out = capture_all %>%
  inner_join(sample_mb %>%
               rename(sample_id_microbiome = sample_id,
                      sample_name_microbiome = sample_name) %>%
               select(-any_of(sample_cols_drop)), by = "capture_id") %>%
  left_join(result_bd %>%
              rename(sample_id_bd = sample_id) %>%
              select(-all_of(sample_cols_drop)), by = "capture_id") %>%  
  left_join(result_amp %>%
              rename(sample_id_amp = sample_id) %>%
              select(-all_of(sample_cols_drop)), by = "capture_id") %>%
  select(
    capture_id,
    taxon_capture,
    sample_id_amp,
    sample_name_amp,
    bd_replicates,	bd_detected,	bd_mean_its1_copies_per_swab,
    total_peptides_ug,
    date,
    site,
    region,
    country, 
    body_temp_c ,
    substrate_temp_c ,
    microhabitat_detailed ,
    water_temp_c,
    p_h, 
    tds_ppm,  
    salinity_ppt,  
    soil_moisture_m3_m3, 
    cloud_cover_percent ,
    dissolved_o2_percent ,
    d1_percent_cover ,
    vegetation_cover_percent, 
    air_pressure_mbar ,
    environmental_latitude,
    environmental_longitude,
    site_elevation_m
  ) %>%
  collect()


#write to file for future use
write.table(capture_sample_out, '~/Documents/GitHub/amphibian_dormancy/ribbitr_bd_amp_muc_envir_data.txt', quote = F, sep='\t', row.names=F)
