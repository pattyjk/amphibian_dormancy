#Pulling data from RIBBiTR database for microbiome samples
library(ribbitrrr)
library(tidyverse)
library(dbplyr)
library(RPostgres) 
library(DBI)

#only if need to reconfigure connection
usethis::edit_r_environ()

#connect to database
dbcon = HopToDB("ribbitr")

# load table "all_tables" from schema "public"
mdt <- tbl(dbcon, Id("public", "all_tables")) %>%
  collect()

# load table "all_columns" from schema "public", filtering to schema "survey_data"
mdc <- tbl(dbcon, Id("public", "all_columns")) %>%
  filter(table_schema == "survey_data") %>%
  collect()

# pointers for all tables of interest for Bd qPCR data
db_bdqpcr = tbl(dbcon, Id("survey_data", "bd_qpcr_results"))
db_sample = tbl(dbcon, Id("survey_data", "sample"))
db_capture = tbl(dbcon, Id("survey_data", "capture"))
db_survey = tbl(dbcon, Id("survey_data", "survey"))
db_visit = tbl(dbcon, Id("survey_data", "visit"))
db_site = tbl(dbcon, Id("survey_data", "site"))
db_region = tbl(dbcon, Id("survey_data", "region"))
db_country = tbl(dbcon, Id("survey_data", "country"))

# left join supporting tables
data_bd_capture = db_bdqpcr %>%
  inner_join(db_sample, by = "sample_id") %>%
  inner_join(db_capture, by = "capture_id") %>%
  left_join(db_survey, by = "survey_id") %>%
  left_join(db_visit, by = "visit_id") %>%
  left_join(db_site, by = "site_id") %>%
  left_join(db_region, by = "region_id") %>%
  left_join(db_country, by = "country_id")

# see what columns are available
colnames(data_bd_capture)

#pull data for use
sql_render(data_bd_capture)
data_bd_capture<-data_bd_capture %>%
  collect()

#write to file for future use
write.table(data_bd_capture, './ribbitr_bd_data.txt', quote = F, sep='\t', row.names=F)
