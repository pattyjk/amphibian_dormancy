# =============================================================================
# RIBBiTR Amphibian Skin Microbiome Analysis Pipeline  v4.1
# Response variables: Total_cells_ml, Per_active, ViableNoSporesTbClCorrected
#
# Changes from v4 — peptides added to core model:
#
#   Peptide defence compounds are biologically central to amphibian skin
#   immunity and microbiome composition. Analysis shows:
#     - Only 11% of peptide variance is between locations (individual signal)
#     - r = 0.35 with total cell density; near-zero with % active
#     - No collinearity issues with remaining core predictors
#     - Missing data is not random: Pennsylvania has ~8% coverage (most
#       species never measured), while Brazil (71%), California (46%), and
#       Panama (54%) have adequate coverage.
#
#   STRATEGY — two parallel core model sets:
#
#   (A) FULL DATASET model (no peptides):
#       Location + temp_seasonality + bd_detected + body_temp_c + (1|site)
#       All 4 locations included. Used for Location-level inference and
#       comparisons that require Pennsylvania.
#
#   (B) PEPTIDE model (Brazil + California + Panama only):
#       Location + temp_seasonality + bd_detected + body_temp_c +
#       log_peptides + (1|site)
#       Pennsylvania excluded explicitly (not silently dropped) due to
#       near-complete absence of peptide measurements. RE reduces to
#       (1|site) within the 3-location subset.
#       Used for all inference about peptide effects.
#
#   Both model sets produce coefficient plots and summary tables.
#   The peptide model is not a sensitivity model — it is a primary model
#   for the subset of locations where peptide data were collected.
#
#   RETAINED from v4:
#     - Species RE dropped (perfectly nested in Location)
#     - temp_seasonality as sole bioclim predictor
#     - snout_vent_length dropped (r=0.88 with body_mass_g)
#     - AIC tests for Location × temp_seasonality and Location × bd_detected
#     - body_mass_g, log_bd_copies in supplementary sensitivity models
#     - Viable cells outlier investigation
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
required_pkgs <- c(
  "tidyverse", "lme4", "lmerTest",
  "MuMIn", "broom.mixed", "performance",
  "ggplot2", "ggeffects", "patchwork",
  "ranger", "ggrepel"
)
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. LOAD & CLEAN DATA
# -----------------------------------------------------------------------------
dat_raw <- read.delim("ribbitr_dormancy_data.txt", sep = "\t",
                      na.strings = c("NA", "N/A", "", "na"),
                      stringsAsFactors = FALSE) %>%
  rename(Total_bacteria_blue = `Total.bacteria_blue`,
         Dilution_factor     = `Dilution.factor`)

num_cols <- c("Total_cells_ml", "Per_active", "ViableNoSporesTbClCorrected",
              "bd_mean_its1_copies_per_swab", "total_peptides_ug",
              "body_temp_c", "body_mass_g", "snout_vent_length",
              "annual_mean_temp", "temp_seasonality", "precip_driest_month",
              "annual_precip", "precip_seasonality", "lat", "long")

dat_raw <- dat_raw %>%
  mutate(
    across(all_of(num_cols), as.numeric),
    bd_detected = case_when(
      bd_detected %in% c("TRUE",  "true",  "1") ~ "Positive",
      bd_detected %in% c("FALSE", "false", "0") ~ "Negative",
      TRUE ~ NA_character_
    ),
    bd_detected = factor(bd_detected, levels = c("Negative", "Positive")),
    Species  = factor(Species),
    Location = factor(Location),
    site     = factor(site)
  ) %>%
  mutate(
    ViableNoSporesTbClCorrected = if_else(
      ViableNoSporesTbClCorrected < 0, NA_real_,
      ViableNoSporesTbClCorrected),
    log_Total_cells_ml = log1p(Total_cells_ml),
    logit_Per_active   = log(Per_active / (1 - Per_active)),
    log_Viable         = log1p(ViableNoSporesTbClCorrected),
    log_bd_copies      = log1p(bd_mean_its1_copies_per_swab),
    log_peptides       = log1p(total_peptides_ug)
  )

dat <- dat_raw %>%
  filter(!is.na(log_Total_cells_ml) | !is.na(logit_Per_active) |
           !is.na(log_Viable))

cat("Total rows:", nrow(dat), "\n")
cat("Location counts:\n"); print(table(dat$Location))
cat("\nSpecies × Location (confirming perfect nesting):\n")
print(table(dat$Species, dat$Location))

# -----------------------------------------------------------------------------
# 2. PREDICTOR STRUCTURE REPORT
# -----------------------------------------------------------------------------
cat("\n=== Env variance: within vs between Location ===\n")
env_check <- c("annual_mean_temp", "temp_seasonality",
               "precip_driest_month", "annual_precip", "precip_seasonality")
for (v in env_check) {
  tot <- var(dat[[v]], na.rm = TRUE)
  win <- mean(tapply(dat[[v]], dat$Location, var, na.rm = TRUE), na.rm = TRUE)
  cat(sprintf("  %-28s  %%var between locations = %3.0f%%\n",
              v, 100 * (1 - win / tot)))
}

cat("\n=== Peptide availability by Location ===\n")
pep_avail <- dat %>%
  group_by(Location) %>%
  summarise(n_total    = n(),
            n_pep      = sum(!is.na(log_peptides)),
            pct_pep    = round(100 * mean(!is.na(log_peptides)), 1),
            .groups    = "drop")
print(pep_avail)

pep_within_var <- dat %>%
  group_by(Location) %>%
  summarise(v = var(log_peptides, na.rm = TRUE), .groups = "drop") %>%
  pull(v) %>% mean(na.rm = TRUE)
pep_total_var <- var(dat$log_peptides, na.rm = TRUE)
cat(sprintf("\n  Peptide %% variance between locations: %.0f%%\n",
            100 * (1 - pep_within_var / pep_total_var)))
cat("  → Peptides are primarily an individual-level signal (within-location)\n")

cat("\n=== Missing data: all predictors ===\n")
for (v in c("temp_seasonality", "precip_driest_month", "bd_detected",
            "body_temp_c", "log_peptides", "body_mass_g", "log_bd_copies")) {
  n_miss <- sum(is.na(dat[[v]]))
  cat(sprintf("  %-22s  missing = %3d  (%4.1f%%)\n",
              v, n_miss, 100 * n_miss / nrow(dat)))
}

# -----------------------------------------------------------------------------
# 3. PREDICTOR SETS
# -----------------------------------------------------------------------------

# (A) FULL DATASET — no peptides (all 4 locations)
fe_full      <- c("Location", "temp_seasonality", "bd_detected", "body_temp_c")
fe_full_env  <- c("Location", "temp_seasonality")
fe_full_host <- c("Location", "bd_detected", "body_temp_c")

# Interaction candidates for full-dataset models (AIC-tested)
fe_full_temp_int <- c(fe_full, "Location:temp_seasonality")
fe_full_bd_int   <- c(fe_full, "Location:bd_detected")

# (B) PEPTIDE MODEL — Brazil + California + Panama only
# Pennsylvania excluded: 92% of PA observations have no peptide data,
# concentrated in 8 species never measured for peptides.
locs_with_pep <- c("Brazil", "California", "Panama")
dat_pep <- dat %>% filter(Location %in% locs_with_pep) %>%
  mutate(Location = droplevels(Location))

fe_pep      <- c("Location", "temp_seasonality", "bd_detected",
                 "body_temp_c", "log_peptides")
fe_pep_env  <- c("Location", "temp_seasonality")
fe_pep_host <- c("Location", "bd_detected", "body_temp_c", "log_peptides")

# Interaction candidate for peptide models
fe_pep_bd_int <- c(fe_pep, "Location:bd_detected")

cat("\n=== Peptide model dataset (Brazil + California + Panama) ===\n")
cat("Rows:", nrow(dat_pep), "\n")
print(table(dat_pep$Location))

# Supplementary sensitivity predictors (remain on full dataset)
fe_sens_mass     <- c(fe_full, "body_mass_g")
fe_sens_bdcopies <- c(fe_full, "log_bd_copies")
fe_sens_full     <- c(fe_full, "body_mass_g", "log_bd_copies")

# RE: site only
re_str <- "(1 | site)"

response_labels <- c(
  "log_Total_cells_ml" = "log(Total cells/mL)",
  "logit_Per_active"   = "logit(% Active)",
  "log_Viable"         = "log(Viable cells, no spores)"
)

# -----------------------------------------------------------------------------
# 4. FIT FUNCTION
# -----------------------------------------------------------------------------
fit_lmm <- function(response, fe_vec, dat_in, re_str, label = NULL) {
  lbl <- if (is.null(label)) response else label
  f   <- as.formula(paste(response, "~",
                          paste(fe_vec, collapse = " + "), "+", re_str))
  fe_cols <- unique(unlist(strsplit(fe_vec, ":")))
  needed  <- unique(c(response, fe_cols,
                      gsub("[^A-Za-z_]", " ",
                           gsub("\\(|\\)", "", re_str)) %>%
                        strsplit(" ") %>% unlist() %>%
                        .[. != "" & . != "1"]))
  needed <- needed[needed %in% names(dat_in)]
  cc <- dat_in %>% select(all_of(needed)) %>% drop_na()
  cat(sprintf("  %-58s n = %d\n", lbl, nrow(cc)))
  m <- tryCatch(
    lmer(f, data = cc, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa",
                               optCtrl   = list(maxfun = 2e5))),
    error = function(e) { warning(lbl, ": ", e$message); NULL }
  )
  if (!is.null(m)) m <- update(m, REML = TRUE)
  list(model = m, data = cc, response = response, label = lbl, n = nrow(cc))
}

aic_select <- function(obj_a, obj_b, label_a = "main", label_b = "interaction") {
  if (is.null(obj_a$model) || is.null(obj_b$model)) return(obj_a)
  ma <- update(obj_a$model, REML = FALSE)
  mb <- update(obj_b$model, REML = FALSE)
  d  <- AIC(ma, mb)
  delta <- d[2, "AIC"] - d[1, "AIC"]
  cat(sprintf("    AIC %s = %.1f  |  %s = %.1f  |  delta = %.2f  →  %s\n",
              label_a, d[1,"AIC"], label_b, d[2,"AIC"], delta,
              if (delta < -2) label_b else label_a))
  if (delta < -2) obj_b else obj_a
}

# -----------------------------------------------------------------------------
# 5. MODEL FITTING — (A) FULL DATASET (all 4 locations, no peptides)
# -----------------------------------------------------------------------------
cat("\n========== MODEL SET A: Full dataset (all locations) ==========\n")

cat("\n--- log(Total cells/mL) ---\n")
mA_cells_core <- fit_lmm("log_Total_cells_ml", fe_full,      dat, re_str, "cells core")
mA_cells_env  <- fit_lmm("log_Total_cells_ml", fe_full_env,  dat, re_str, "cells env-only")
mA_cells_host <- fit_lmm("log_Total_cells_ml", fe_full_host, dat, re_str, "cells host-only")
mA_cells_int  <- fit_lmm("log_Total_cells_ml", fe_full_temp_int, dat, re_str,
                         "cells [+Location×temp_seas]")
cat("  AIC: core vs Location×temp_seasonality:\n")
bestA_cells <- aic_select(mA_cells_core, mA_cells_int)

cat("\n--- logit(% Active) ---\n")
mA_active_core <- fit_lmm("logit_Per_active", fe_full,      dat, re_str, "active core")
mA_active_env  <- fit_lmm("logit_Per_active", fe_full_env,  dat, re_str, "active env-only")
mA_active_host <- fit_lmm("logit_Per_active", fe_full_host, dat, re_str, "active host-only")
mA_active_int  <- fit_lmm("logit_Per_active", fe_full_bd_int, dat, re_str,
                          "active [+Location×bd_detected]")
cat("  AIC: core vs Location×bd_detected:\n")
bestA_active <- aic_select(mA_active_core, mA_active_int)

cat("\n--- log(Viable cells, no spores) ---\n")
mA_viable_core <- fit_lmm("log_Viable", fe_full,      dat, re_str, "viable core")
mA_viable_env  <- fit_lmm("log_Viable", fe_full_env,  dat, re_str, "viable env-only")
mA_viable_host <- fit_lmm("log_Viable", fe_full_host, dat, re_str, "viable host-only")

model_list_A <- list(
  log_Total_cells_ml = list(global = bestA_cells$model,  env = mA_cells_env$model,
                            host   = mA_cells_host$model, data = bestA_cells$data),
  logit_Per_active   = list(global = bestA_active$model, env = mA_active_env$model,
                            host   = mA_active_host$model, data = bestA_active$data),
  log_Viable         = list(global = mA_viable_core$model, env = mA_viable_env$model,
                            host   = mA_viable_host$model, data = mA_viable_core$data)
)

# -----------------------------------------------------------------------------
# 6. MODEL FITTING — (B) PEPTIDE MODELS (Brazil + California + Panama)
# -----------------------------------------------------------------------------
cat("\n========== MODEL SET B: Peptide models (3 locations) ==========\n")
cat("  Pennsylvania excluded: peptide data absent for 8/9 species\n")
cat("  n available:", nrow(dat_pep), "rows\n\n")

cat("--- log(Total cells/mL) ---\n")
mB_cells_core <- fit_lmm("log_Total_cells_ml", fe_pep,      dat_pep, re_str,
                         "cells+pep core")
mB_cells_env  <- fit_lmm("log_Total_cells_ml", fe_pep_env,  dat_pep, re_str,
                         "cells+pep env-only")
mB_cells_host <- fit_lmm("log_Total_cells_ml", fe_pep_host, dat_pep, re_str,
                         "cells+pep host-only")

cat("\n--- logit(% Active) ---\n")
mB_active_core <- fit_lmm("logit_Per_active", fe_pep,      dat_pep, re_str,
                          "active+pep core")
mB_active_env  <- fit_lmm("logit_Per_active", fe_pep_env,  dat_pep, re_str,
                          "active+pep env-only")
mB_active_host <- fit_lmm("logit_Per_active", fe_pep_host, dat_pep, re_str,
                          "active+pep host-only")
mB_active_int  <- fit_lmm("logit_Per_active", fe_pep_bd_int, dat_pep, re_str,
                          "active+pep [+Location×bd_detected]")
cat("  AIC: core vs Location×bd_detected:\n")
bestB_active <- aic_select(mB_active_core, mB_active_int)

cat("\n--- log(Viable cells, no spores) ---\n")
mB_viable_core <- fit_lmm("log_Viable", fe_pep,      dat_pep, re_str,
                          "viable+pep core")
mB_viable_env  <- fit_lmm("log_Viable", fe_pep_env,  dat_pep, re_str,
                          "viable+pep env-only")
mB_viable_host <- fit_lmm("log_Viable", fe_pep_host, dat_pep, re_str,
                          "viable+pep host-only")

model_list_B <- list(
  log_Total_cells_ml = list(global = mB_cells_core$model,   env = mB_cells_env$model,
                            host   = mB_cells_host$model,   data = mB_cells_core$data),
  logit_Per_active   = list(global = bestB_active$model,     env = mB_active_env$model,
                            host   = mB_active_host$model,  data = bestB_active$data),
  log_Viable         = list(global = mB_viable_core$model,   env = mB_viable_env$model,
                            host   = mB_viable_host$model,  data = mB_viable_core$data)
)

# Also keep a single reference for viable cells outlier section
m_viable_core <- mA_viable_core

# -----------------------------------------------------------------------------
# 6b. SUPPLEMENTARY SENSITIVITY MODELS (full dataset, optional predictors)
# -----------------------------------------------------------------------------
cat("\n--- Supplementary sensitivity (full dataset) ---\n")
sens_models <- lapply(
  list(
    list(rsp = "log_Total_cells_ml", fe = fe_sens_mass,     lbl = "cells + body_mass"),
    list(rsp = "log_Total_cells_ml", fe = fe_sens_bdcopies, lbl = "cells + bd_copies"),
    list(rsp = "log_Total_cells_ml", fe = fe_sens_full,     lbl = "cells + body_mass + bd_copies"),
    list(rsp = "logit_Per_active",   fe = fe_sens_mass,     lbl = "active + body_mass"),
    list(rsp = "logit_Per_active",   fe = fe_sens_bdcopies, lbl = "active + bd_copies"),
    list(rsp = "logit_Per_active",   fe = fe_sens_full,     lbl = "active + body_mass + bd_copies"),
    list(rsp = "log_Viable",         fe = fe_sens_full,     lbl = "viable + body_mass + bd_copies")
  ),
  function(x) fit_lmm(x$rsp, x$fe, dat, re_str, x$lbl)
)

# -----------------------------------------------------------------------------
# 7. MODEL SUMMARY TABLES
# -----------------------------------------------------------------------------
nice_labels <- c(
  "LocationCalifornia"                        = "Location [California]",
  "LocationPanama"                            = "Location [Panama]",
  "LocationPennsylvania"                      = "Location [Pennsylvania]",
  "temp_seasonality"                          = "Temp seasonality",
  "bd_detectedPositive"                       = "Bd detected [+]",
  "body_temp_c"                               = "Body temperature (°C)",
  "log_peptides"                              = "log(Peptides µg)",
  "log_bd_copies"                             = "log(Bd ITS1 copies)",
  "body_mass_g"                               = "Body mass (g)",
  "LocationCalifornia:temp_seasonality"       = "California × Temp seasonality",
  "LocationPanama:temp_seasonality"           = "Panama × Temp seasonality",
  "LocationPennsylvania:temp_seasonality"     = "Pennsylvania × Temp seasonality",
  "LocationCalifornia:bd_detectedPositive"    = "California × Bd detected [+]",
  "LocationPanama:bd_detectedPositive"        = "Panama × Bd detected [+]",
  "LocationPennsylvania:bd_detectedPositive"  = "Pennsylvania × Bd detected [+]"
)

extract_table <- function(model, label) {
  if (is.null(model)) return(NULL)
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      response = label,
      across(where(is.numeric), ~ round(.x, 4)),
      sig = case_when(
        p.value < 0.001 ~ "***", p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",   p.value < 0.1   ~ ".", TRUE ~ "")
    ) %>%
    select(response, term, estimate, std.error, conf.low, conf.high,
           statistic, p.value, sig)
}

make_tables <- function(ml, suffix) {
  lapply(names(ml), function(rsp) {
    bind_rows(
      extract_table(ml[[rsp]]$global, paste0(rsp, " [", suffix, " global]")),
      extract_table(ml[[rsp]]$env,    paste0(rsp, " [", suffix, " env-only]")),
      extract_table(ml[[rsp]]$host,   paste0(rsp, " [", suffix, " host-only]"))
    )
  }) %>% bind_rows()
}

tables_A <- make_tables(model_list_A, "full")
tables_B <- make_tables(model_list_B, "peptide")

sens_tables <- lapply(sens_models, function(obj) {
  extract_table(obj$model, obj$label)
}) %>% bind_rows()

all_tables <- bind_rows(tables_A, tables_B, sens_tables)
write.csv(all_tables, "model_summary_tables_v4.1.csv", row.names = FALSE)
cat("\nAll tables saved: model_summary_tables_v4.1.csv\n")

# Separate files for convenience
write.csv(tables_A, "model_summary_fulldat_v4.1.csv",    row.names = FALSE)
write.csv(tables_B, "model_summary_peptides_v4.1.csv",   row.names = FALSE)
write.csv(sens_tables, "model_summary_sensitivity_v4.1.csv", row.names = FALSE)
cat("Individual table files saved.\n")

# R² for both model sets
r2_tab <- lapply(list(A = model_list_A, B = model_list_B), function(ml) {
  lapply(names(ml), function(rsp) {
    m <- ml[[rsp]]$global
    if (is.null(m)) return(NULL)
    r2 <- MuMIn::r.squaredGLMM(m)
    data.frame(model_set = if (identical(ml, model_list_A)) "full" else "peptide",
               response  = response_labels[rsp],
               R2m = round(r2[1, "R2m"], 3),
               R2c = round(r2[1, "R2c"], 3))
  }) %>% bind_rows()
}) %>% bind_rows()
cat("\nMarginal / Conditional R²:\n"); print(r2_tab)

# -----------------------------------------------------------------------------
# 8. RANDOM FOREST — both model sets
# -----------------------------------------------------------------------------
fit_rf <- function(response, fe_vec, dat_in) {
  preds <- unique(c(unlist(strsplit(fe_vec, ":")), "site"))
  preds <- preds[preds %in% names(dat_in)]
  cc    <- dat_in %>%
    select(all_of(c(response, preds))) %>% drop_na() %>%
    mutate(across(where(is.factor), as.integer))
  tryCatch(
    ranger(as.formula(paste(response, "~ .")),
           data = cc, num.trees = 1000,
           importance = "permutation", seed = 42),
    error = function(e) NULL
  )
}

build_rf_imp <- function(fe_vec, dat_in, set_label) {
  bind_rows(lapply(names(response_labels), function(rsp) {
    rf <- fit_rf(rsp, fe_vec, dat_in)
    if (is.null(rf)) return(NULL)
    data.frame(predictor  = names(importance(rf)),
               importance = importance(rf),
               resp_label = response_labels[rsp],
               model_set  = set_label)
  })) %>%
    mutate(driver_type = if_else(
      str_detect(predictor, paste(c("bd","body","peptide"), collapse = "|")),
      "Host", "Environment"))
}

rf_imp_A <- build_rf_imp(fe_full, dat,     "Full dataset (4 locations)")
rf_imp_B <- build_rf_imp(fe_pep,  dat_pep, "Peptide model (3 locations)")
rf_imp   <- bind_rows(rf_imp_A, rf_imp_B)

rf_plot <- rf_imp %>%
  group_by(model_set, resp_label) %>%
  slice_max(importance, n = 8) %>%
  ggplot(aes(importance, reorder(predictor, importance), fill = driver_type)) +
  geom_col(colour = "white", linewidth = 0.3) +
  facet_grid(model_set ~ resp_label, scales = "free") +
  scale_fill_manual(values = c("Host" = "#E07B54", "Environment" = "#4A90A4"),
                    name = "Driver type") +
  labs(title = "Random Forest: permutation variable importance",
       x = "Mean decrease in MSE", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 8))

ggsave("rf_variable_importance_v4.1.png", rf_plot,
       width = 13, height = 8, dpi = 200)
cat("RF importance plot saved: rf_variable_importance_v4.1.png\n")

# -----------------------------------------------------------------------------
# 9. COEFFICIENT PLOTS — both model sets
# -----------------------------------------------------------------------------
build_coef_data <- function(ml, set_label) {
  lapply(names(ml), function(rsp) {
    m <- ml[[rsp]]$global
    if (is.null(m)) return(NULL)
    broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        response       = rsp,
        resp_label     = response_labels[rsp],
        model_set      = set_label,
        term_clean     = if_else(term %in% names(nice_labels),
                                 nice_labels[term], term),
        driver_type    = if_else(str_detect(term, "bd|body|peptide"),
                                 "Host", "Environment"),
        is_interaction = str_detect(term, ":")
      )
  }) %>% bind_rows()
}

coef_A <- build_coef_data(model_list_A, "Full dataset (4 locations)")
coef_B <- build_coef_data(model_list_B, "Peptide model (3 locations)")

make_coef_panel <- function(coef_df, set_label, subtitle_extra = "") {
  plots <- lapply(unique(coef_df$response), function(rsp) {
    d <- coef_df %>% filter(response == rsp)
    if (nrow(d) == 0) return(NULL)
    ggplot(d, aes(estimate, reorder(term_clean, estimate),
                  colour = driver_type, shape = is_interaction)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                     height = 0.25, linewidth = 0.7) +
      geom_point(size = 3) +
      scale_colour_manual(
        values = c("Host" = "#E07B54", "Environment" = "#4A90A4"),
        name   = "Driver type") +
      scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 17),
        labels = c("Main effect", "Interaction"), name = NULL) +
      labs(title = unique(d$resp_label),
           x = "Coefficient estimate (± 95% CI)", y = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold"))
  })
  plots <- plots[!sapply(plots, is.null)]
  wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title    = paste("Fixed-effect coefficients:", set_label),
      subtitle = paste0("RE: (1|site)", subtitle_extra),
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )
}

p_coef_A <- make_coef_panel(coef_A, "Full dataset (all 4 locations)",
                            " | No peptides; Pennsylvania included")
p_coef_B <- make_coef_panel(coef_B, "Peptide models (Brazil + California + Panama)",
                            " | Pennsylvania excluded (insufficient peptide data)")

ggsave("coef_plot_fulldat_v4.1.png",   p_coef_A, width = 9, height = 13, dpi = 200)
ggsave("coef_plot_peptides_v4.1.png",  p_coef_B, width = 9, height = 13, dpi = 200)
cat("Coefficient plots saved: coef_plot_fulldat_v4.1.png, coef_plot_peptides_v4.1.png\n")

# -----------------------------------------------------------------------------
# 10. MODEL DIAGNOSTICS — both model sets
# -----------------------------------------------------------------------------
loc_cols <- c("Brazil"       = "#E07B54", "California"   = "#4A90A4",
              "Panama"       = "#8E6BBF", "Pennsylvania" = "#5B9E6E")

make_diag_panel <- function(ml, title_str) {
  diag_plots <- lapply(names(ml), function(rsp) {
    m <- ml[[rsp]]$global
    if (is.null(m)) return(NULL)
    idx    <- as.integer(names(fitted(m)))
    res_df <- data.frame(
      fitted   = as.numeric(fitted(m)),
      residual = resid(m, type = "pearson"),
      Location = ml[[rsp]]$data$Location[idx]
    )
    p_rf <- ggplot(res_df, aes(fitted, residual)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(colour = Location), alpha = 0.45, size = 1.4) +
      geom_smooth(se = FALSE, colour = "red", linewidth = 0.8) +
      scale_colour_manual(values = loc_cols) +
      labs(title    = response_labels[rsp], subtitle = "Colour = Location",
           x = "Fitted", y = "Pearson residual") +
      theme_bw(base_size = 10) + theme(legend.position = "bottom")
    p_qq <- ggplot(res_df, aes(sample = residual)) +
      stat_qq(aes(colour = Location), alpha = 0.45) +
      stat_qq_line(colour = "red") +
      scale_colour_manual(values = loc_cols) +
      labs(title = response_labels[rsp],
           x = "Theoretical quantiles", y = "Sample quantiles") +
      theme_bw(base_size = 10) + theme(legend.position = "bottom")
    list(resfit = p_rf, qq = p_qq)
  })
  wrap_plots(
    c(lapply(diag_plots, `[[`, "resfit"),
      lapply(diag_plots, `[[`, "qq")),
    ncol = 3
  ) + plot_annotation(
    title    = title_str,
    subtitle = "Flat smoother + interspersed location colours = well-specified model",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 9,  colour = "grey40")))
}

diag_A <- make_diag_panel(model_list_A, "Diagnostics: Full dataset models (v4.1)")
diag_B <- make_diag_panel(model_list_B, "Diagnostics: Peptide models — 3 locations (v4.1)")

ggsave("model_diagnostics_fulldat_v4.1.png",  diag_A, width = 14, height = 9, dpi = 200)
ggsave("model_diagnostics_peptides_v4.1.png", diag_B, width = 14, height = 9, dpi = 200)
cat("Diagnostic plots saved.\n")

# -----------------------------------------------------------------------------
# 11. VIABLE CELLS — OUTLIER INVESTIGATION (full dataset model)
# -----------------------------------------------------------------------------
m_viable        <- mA_viable_core$model
dat_viable_used <- mA_viable_core$data

if (!is.null(m_viable)) {
  idx      <- as.integer(names(fitted(m_viable)))
  dat_used <- dat_viable_used[idx, ]
  std_r    <- resid(m_viable, type = "pearson") /
    sd(resid(m_viable, type = "pearson"))
  
  outlier_df <- dat_used %>%
    mutate(std_resid  = std_r,
           fitted_val = as.numeric(fitted(m_viable))) %>%
    filter(abs(std_resid) > 3) %>%
    select(RIBBITR_ID, Species, Location, site,
           ViableNoSporesTbClCorrected, log_Viable, std_resid) %>%
    arrange(std_resid)
  
  cat(sprintf("\n--- Viable cells outliers (|std resid| > 3): %d / %d ---\n",
              nrow(outlier_df), nrow(dat_used)))
  print(outlier_df)
  write.csv(outlier_df, "viable_cells_outliers_v4.1.csv", row.names = FALSE)
  
  # Sensitivity: remove outliers
  if (nrow(outlier_df) > 0) {
    dat_noout <- dat %>% filter(!RIBBITR_ID %in% outlier_df$RIBBITR_ID)
    cat("\nSensitivity model (outliers removed):\n")
    m_viable_noout <- fit_lmm("log_Viable", fe_full, dat_noout,
                              re_str, "viable [outliers removed]")
    sens_comp <- bind_rows(
      extract_table(m_viable,               "viable [original]"),
      extract_table(m_viable_noout$model,   "viable [outliers removed]")
    ) %>% filter(term != "(Intercept)")
    write.csv(sens_comp, "viable_sensitivity_v4.1.csv", row.names = FALSE)
    cat("Viable sensitivity comparison saved: viable_sensitivity_v4.1.csv\n")
  }
  
  # Outlier plot
  full_df <- dat_used %>%
    mutate(std_resid  = std_r,
           fitted_val = as.numeric(fitted(m_viable)),
           is_outlier = abs(std_r) > 3)
  
  p_out <- ggplot(full_df, aes(fitted_val, std_resid)) +
    geom_hline(yintercept = c(-3, 0, 3),
               linetype   = c("dashed","solid","dashed"),
               colour     = c("red","grey40","red")) +
    geom_point(aes(colour = Location, shape = is_outlier),
               alpha = 0.55, size = 1.8) +
    ggrepel::geom_text_repel(
      data   = full_df %>% filter(is_outlier),
      aes(label = RIBBITR_ID), size = 2.5, colour = "red",
      max.overlaps = 20) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_colour_manual(
      values = c("Brazil"       = "#E07B54", "California"   = "#4A90A4",
                 "Panama"       = "#8E6BBF", "Pennsylvania" = "#5B9E6E")) +
    labs(title    = "log(Viable cells) — outlier identification (v4.1)",
         subtitle = "Triangles and labels = |std resid| > 3",
         x = "Fitted values", y = "Standardised Pearson residual") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave("viable_cells_outlier_plot_v4.1.png", p_out,
         width = 8, height = 5, dpi = 200)
  cat("Outlier plot saved: viable_cells_outlier_plot_v4.1.png\n")
}

# -----------------------------------------------------------------------------
# 12. SESSION INFO
# -----------------------------------------------------------------------------
cat("\n--- Session Info ---\n")
sessionInfo()