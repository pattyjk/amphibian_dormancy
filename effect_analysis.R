#look at corerlations
# Load necessary libraries
# install.packages(c("corrplot", "dplyr", "reshape2", "ggplot2"))
library(corrplot)
library(dplyr)
library(reshape2)
library(ggplot2)

# 2. Re-import specifying multiple types of "Not Available" strings
df <- read.delim("~/Documents/GitHub/amphibian_dormancy/ribbitr_dorm_pruned.txt", 
                 sep="\t", 
                 header=TRUE, 
                 na.strings = c("NA", "N/A", "", " "), 
                 stringsAsFactors = FALSE)
#remove tube id cols
df<-df[,-c(11:13)]

# 3. List the columns you want to analyze (Biological + Environmental)
# We use names(df) to grab all columns from snout_vent_length to the end
metric_names <- names(df)[9:ncol(df)]

# 4. Force these columns to be numeric (this fixes the "text" issue)
numeric_df <- df %>%
  mutate(across(all_of(metric_names), as.numeric)) %>%
  select(where(is.numeric))

# 5. Check the structure again - you should now see ~40 variables
str(numeric_df)

# 6. Generate the Correlation Plot
# 'pairwise.complete.obs' is critical so the plot doesn't turn out blank/white
cor_matrix <- cor(numeric_df, use = "pairwise.complete.obs")

# Create the plot
# If you are in RStudio, this will appear in the 'Plots' pane
corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.5,
         mar = c(0,0,1,0),
         title = "Amphibian Microbial Activity Correlations")

###############
library(ggplot2)
library(dplyr)

# 1. Regional and Species comparison
# Filter for top 10 most common species to keep the plot clean
top_species <- df %>% count(Species, sort=TRUE) %>% top_n(10) %>% pull(Species)
df_top <- df %>% filter(Species %in% top_species)

ggplot(df_top, aes(x=Location, y=Per_active, fill=Species)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title="Microbial Activity by Location and Species")

# ANOVA test
anova_res <- aov(Per_active ~ Location + Species, data=df)
summary(anova_res)
TukeyHSD(anova_res)

# 2. Environmental Stress vs Dormancy
# Regression Plot for Seasonality
ggplot(df, aes(x=annual_mean_temp, y=Per_dormant)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", color="blue") +
  labs(title="Seasonality vs Microbial Dormancy",
       x="Annual Mean Temp", y="Percent Dormant")

# Linear Model
model_dormancy <- lm(Per_dormant ~ mean_temp_coldest_quarter + precip_driest_month, data=df)
summary(model_dormancy)
#########
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)

# 1. Clean and Prepare
# Assuming 'df' is your loaded data
df_clean <- df %>%
  mutate(across(c(Per_active, PerActive_Ratio, Total_cells_ml, TotalNoSporesTbClCorrected,
                  snout_vent_length, body_mass_g, bd_mean_its1_copies_per_swab, 
                  lat, long, annual_mean_temp, temp_seasonality, 
                  annual_precip, precip_seasonality), as.numeric)) %>%
  mutate(bd_presence = as.numeric(as.logical(bd_detected)))

# 2. Function to run standardized models
run_std_model <- function(target_var, data) {
  # Standardize all variables for that specific model
  data_subset <- data %>%
    select(all_of(target_var), snout_vent_length, body_mass_g, bd_presence, 
           lat, long, annual_mean_temp, temp_seasonality, 
           annual_precip, precip_seasonality, bd_mean_its1_copies_per_swab) %>%
    drop_na() %>%
    mutate(across(everything(), scale))
  
  form <- as.formula(paste(target_var, "~ ."))
  lm(form, data = data_subset) %>% 
    tidy(conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(target = target_var)
}

# 3. Combine results for all 4 targets
targets <- c("Per_active", "PerActive_Ratio", "Total_cells_ml", "TotalNoSporesTbClCorrected")
all_results <- bind_rows(lapply(targets, run_std_model, data = df_clean))

# 4. Create the Forest Plot
cb_palette <- c(
  "Per_active" = "#E69F00",                # Orange
  "PerActive_Ratio" = "#56B4E9",           # Sky Blue
  "Total_cells_ml" = "#009E73",            # Bluish Green
  "TotalNoSporesTbClCorrected" = "#D55E00" # Vermillion
)

#Create the Forest Plot for all
ggplot(all_results, aes(x = estimate, y = term, color = target)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, position = position_dodge(width = 0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  
  # Apply the color-blind friendly palette and custom labels
  scale_color_manual(
    name = "Microbial Metric",
    values = cb_palette,
    labels = c(
      "Per_active" = "Percent active (Cell staining)",
      "PerActive_Ratio" = "Percent active (Taxonomic ratio)",
      "Total_cells_ml" = "Total microbial load (Cells/mL)",
      "TotalNoSporesTbClCorrected" = "Total endospore count"
    )
  ) +
  
  labs(
    title = "Standardized Drivers of Amphibian Microbial States",
    subtitle = "Standardized Beta Coefficients with 95% Confidence Intervals",
    x = "Effect Size (Standardized Beta)", 
    y = "Predictor Variable"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )
###########
#sep by host/pathogen and environmental
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)

# 1. Prepare the Color-Blind Palette and Labels
cb_palette <- c(
  "Per_active" = "#E69F00",                # Orange
  "PerActive_Ratio" = "#56B4E9",           # Sky Blue
  "Total_cells_ml" = "#009E73",            # Bluish Green
  "TotalNoSporesTbClCorrected" = "#D55E00" # Vermillion
)

target_labels <- c(
  "Per_active" = "Percent active (Cell staining)",
  "PerActive_Ratio" = "Percent active (16S ratio)",
  "Total_cells_ml" = "Total microbial load (Cells/mL)",
  "TotalNoSporesTbClCorrected" = "Total endospore count"
)

# 2. Run Standardized Regression (Run one model with ALL predictors)
# [Assume results are stored in 'all_results' from previous steps]

# 3. Define Predictor Groups
host_pathogen_vars <- c("snout_vent_length", "body_mass_g", "bd_presence", 
                        "sexMale", "sexUnknown", "lifeJuvenile", "lifeSubadult", 'bd_mean_its1_copies_per_swab')
env_vars <- c("lat", "long", "annual_mean_temp", "temp_seasonality", 
              "annual_precip", "precip_seasonality")

predictor_labels <- c(
  "snout_vent_length" = "Body Size (SVL)",
  "body_mass_g" = "Body Mass (g)",
  "bd_presence" = "Pathogen (Bd) Presence",
  "bd_mean_its1_copies_per_swab" = "Infection Intensity (ITS1 abundance)",
  "lat" = "Latitude",
  "long" = "Longitude",
  "annual_mean_temp" = "Annual Mean Temp",
  "temp_seasonality" = "Temp. Varibility (Temp. Seasonality)",
  "annual_precip" = "Annual Precipitation (mm)",
  "precip_seasonality" = "Precipitation Seasonality",
  "sex_Male" = "Sex (Male)",
  "sex_Unknown" = "Sex (Unknown)",
  "life_Juvenile" = "Life Stage (Juvenile)",
  "life_Subadult" = "Life Stage (Subadult)"
)

metric_labels <- c(
  "Per_active" = "Percent active (Cell staining)",
  "PerActive_Ratio" = "Percent active (Taxonomic ratio)",
  "Total_cells_ml" = "Total microbial load (Cells/mL)",
  "TotalNoSporesTbClCorrected" = "Total spore count"
)

# 2. Build the Plot Function
plot_forest_clean <- function(data_subset, title) {
  ggplot(data_subset, aes(x = estimate, y = term, color = target)) +
    geom_point(position = position_dodge(width = 0.7), size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, position = position_dodge(width = 0.7)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    # Apply Color-Blind Palette and Metric Labels
    scale_color_manual(values = cb_palette, labels = metric_labels, name = "Metric") +
    
    # Apply Meaningful Predictor Labels
    scale_y_discrete(labels = predictor_labels) +
    
    theme_minimal() +
    labs(title = title, x = "Standardized Impact (Beta)", y = "Predictor Variable")
}

# 3. Generate the Plots
# Assuming 'all_results' contains the regression output
plot_forest_clean(all_results %>% filter(term %in% host_pathogen_vars), "Host & Pathogen Drivers")
plot_forest_clean(all_results %>% filter(term %in% env_vars), "Environmental Drivers")

###get stats
library(dplyr)
library(tidyr)
library(broom)

# 1. Load and Clean Data
df <- read.delim("~/Documents/GitHub/amphibian_dormancy/ribbitr_dorm_pruned.txt", 
                 sep="\t", 
                 header=TRUE, 
                 na.strings = c("NA", "N/A", "", " "), 
                 stringsAsFactors = FALSE)

# 2. Prepare Variables
# Convert bd_detected to a 0/1 numeric and standardize numeric predictors
df_clean <- df %>%
  mutate(
    bd_presence = as.numeric(as.logical(bd_detected)),
    across(c(snout_vent_length, body_mass_g, lat, long, 
             annual_mean_temp, temp_seasonality, 
             annual_precip, precip_seasonality), scale)
  )

# 3. Define the Function to Run and Tidy Models
get_model_stats <- function(target_var) {
  # Build the formula string
  formula_str <- paste(target_var, "~ snout_vent_length + body_mass_g + bd_presence + 
                       lat + long + annual_mean_temp + temp_seasonality + 
                       annual_precip + precip_seasonality + sex + life_stage")
  
  # Run the model
  # We scale the target variable as well to get true Standardized Betas
  model <- lm(as.formula(formula_str), data = df_clean %>% mutate(!!target_var := scale(get(target_var))))
  
  # Extract stats using broom
  stats_table <- tidy(model, conf.int = TRUE) %>%
    mutate(
      target_metric = target_var,
      significant = ifelse(p.value < 0.05, "*", "")
    )
  
  return(stats_table)
}

# 4. Generate Tables for all 4 Metrics
targets <- c("Per_active", "PerActive_Ratio", "Total_cells_ml", "TotalNoSporesTbClCorrected")
all_stats_table <- bind_rows(lapply(targets, get_model_stats))

# 5. View and Export
# Sort by significance and target
all_stats_table <- all_stats_table %>%
  select(target_metric, term, estimate, std.error, p.value, conf.low, conf.high, significant) %>%
  arrange(target_metric, p.value)

# Print to console
print(all_stats_table)

# Export to CSV 
write.csv(all_stats_table, "~/Documents/GitHub/amphibian_dormancy/microbial_analysis_stats_table.csv", row.names = FALSE)
