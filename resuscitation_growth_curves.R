# =============================================================================
# Growth curves for the 4 resuscitation isolates x the 4 shared additions
#
# Inputs : growth_curves_resuscitation_subset.txt
#            (trimmed raw OD600 file: 4 isolates x 4 treatments = 144 curves,
#             0-72 h every 30 min, 1 replicate/curve; same columns as the full
#             growth_curves.txt, which can also be used by changing growth_file)
#          resucitation_data.txt    (resuscitation experiment; for the comparison)
# Outputs: fig_growth_curves.png, fig_growth_dose_response.png,
#          fig_growth_vs_resuscitation.png,
#          growth_curve_parameters.csv, growth_ic50.csv, growth_vs_resuscitation.csv
#
# Per curve: logistic fit N(t) = K / (1 + ((K - N0)/N0) exp(-r t)) by nlsLM
#   (same model as growthcurver), giving K, r, mu_max = max d ln N / dt,
#   plus empirical AUC (trapezoid on raw OD600). Curves are pure logistic from
#   OD ~0.01 with no lag phase, so lag is not reported.
# Dose response: AUC relative to the 0-concentration curve of the same
#   isolate x treatment; IC50 from a log-logistic fit (drc LL.4, top = 1, bottom = 0)
#   where growth falls below 50 %.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(minpack.lm); library(drc)
})
select <- dplyr::select   # drc loads MASS, which masks select()

growth_file <- "~/Documents/GitHub/amphibian_dormancy/growth_curves_resuscitation.txt"   # or the full "growth_curves.txt"
resus_file  <- "~/Documents/GitHub/amphibian_dormancy/resucitation_data.txt"
outdir      <- "~/Documents/GitHub/amphibian_dormancy/"

isolates <- c("126", "343", "336", "616")
trts <- c(NaturalMix_PBS = "Skin peptide mix", Bufadienolides_EtOH = "Bufadienolides",
          Corticosterone_EtOH = "Corticosterone", ChytridMetabolites = "Bd metabolites")

# ---- 1. Read & subset ---------------------------------------------------------
G <- read.delim(growth_file, check.names = FALSE, colClasses = c(isolate_id = "character")) %>%
  filter(isolate_id %in% isolates, treatment %in% names(trts)) %>%
  mutate(species = trimws(sub(" strain.*$", "", species)),
         host = recode(frog_species, RaCa = "Rana catesbeiana", HyPh = "Hylodes phyllodes"),
         concentration_unit = ifelse(concentration_unit == "relative", "relative", "\u00b5M"))   # file uses Greek mu; plot font needs micro sign
stopifnot(n_distinct(G$isolate_id) == 4, n_distinct(G$treatment) == 4)
n_curves <- nrow(distinct(G, isolate_id, treatment, concentration))
if (n_curves != 144) warning("Expected 144 curves (4 isolates x 4 treatments x all doses); found ", n_curves)
cat("Read", basename(growth_file), ":", nrow(G), "rows,", n_curves, "curves\n")
iso_lab <- G %>% distinct(isolate_id, species, host) %>%
  mutate(isolate_id = factor(isolate_id, levels = isolates)) %>% arrange(isolate_id) %>%
  mutate(lab = paste0(species, " (", isolate_id, ")"))
iso_col <- setNames(c("#0072B2", "#56B4E9", "#D55E00", "#E69F00"), isolates)   # Acinetobacter blues, Paenibacillus oranges

# ---- 2. Logistic fit per curve -------------------------------------------------
fit_logistic <- function(t, y) {
  f <- tryCatch(nlsLM(y ~ K / (1 + ((K - N0) / N0) * exp(-r * t)),
                      start = list(K = max(y), N0 = max(y[1], 1e-4), r = 0.3),
                      lower = c(1e-4, 1e-5, 1e-4), upper = c(5, 0.5, 5),
                      control = nls.lm.control(maxiter = 500)), error = function(e) NULL)
  auc <- sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
  if (is.null(f)) return(tibble(K = NA, N0 = NA, r = NA, mu_max = NA, rmse = NA, auc = auc, od_max = max(y)))
  p  <- coef(f); tt <- seq(0, max(t), 0.01)
  lnN <- log(p["K"] / (1 + ((p["K"] - p["N0"]) / p["N0"]) * exp(-p["r"] * tt)))
  tibble(K = p["K"], N0 = p["N0"], r = p["r"], mu_max = max(diff(lnN) / 0.01),
         rmse = sqrt(mean(resid(f)^2)), auc = auc, od_max = max(y))
}
P <- G %>% arrange(timepoint_hours) %>%
  group_by(isolate_id, species, host, treatment, concentration, concentration_unit,
           file_mu = growth_rate_mu, file_K = carrying_capacity_K) %>%
  group_modify(~ fit_logistic(.x$timepoint_hours, .x$od600)) %>% ungroup() %>%
  group_by(isolate_id, treatment) %>%
  mutate(rel_auc = auc / auc[concentration == 0], rel_r = r / r[concentration == 0], rel_K = K / K[concentration == 0]) %>%
  ungroup()
cat("Curves fitted:", nrow(P), " failed:", sum(is.na(P$K)), "\n")
cat("Agreement with file parameters: r vs growth_rate_mu", round(cor(P$r, P$file_mu), 3),
    "; K vs carrying_capacity_K", round(cor(P$K, P$file_K), 3), "\n")
write.csv(P, file.path(outdir, "growth_curve_parameters.csv"), row.names = FALSE)

# ---- 3. IC50 for treatments that reach < 50 % growth ----------------------------
inhib <- P %>% group_by(isolate_id, treatment) %>% filter(min(rel_auc) < 0.5) %>% ungroup()
IC <- inhib %>% group_by(isolate_id, treatment) %>% group_modify(function(x, k) {
  m <- drm(rel_auc ~ concentration, data = x, fct = LL.4(fixed = c(NA, 0, 1, NA)))
  e <- ED(m, 50, interval = "delta", display = FALSE)
  tibble(ic50 = e[1], ic50_lo = e[3], ic50_hi = e[4], hill = unname(coef(m)[1]))
}) %>% ungroup()
write.csv(IC, file.path(outdir, "growth_ic50.csv"), row.names = FALSE)
print(IC %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ---- 4. Shared styling -----------------------------------------------------------
th <- theme_classic(base_size = 9) +
  theme(axis.text = element_text(colour = "black"), plot.title = element_text(size = 9.5),
        plot.title.position = "plot", strip.background = element_blank(),
        strip.text = element_text(size = 8.5), plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30"))
P <- P %>% mutate(trt = factor(trts[treatment], levels = trts),
                  iso = factor(isolate_id, levels = isolates, labels = iso_lab$lab))
G <- G %>% mutate(trt = factor(trts[treatment], levels = trts),
                  iso = factor(isolate_id, levels = isolates, labels = iso_lab$lab))
unit_lab <- G %>% distinct(trt, concentration_unit) %>% arrange(trt)

# ---- 5. Figure 1: raw growth curves (isolate x treatment) -------------------------
# colour = dose rank within treatment (doses differ between treatments); 0 = black
G <- G %>% group_by(treatment) %>% mutate(dose_rank = dense_rank(concentration) - 1,
                                          dose_frac = dose_rank / max(dose_rank)) %>% ungroup()
dose_keys <- G %>% distinct(trt, concentration, concentration_unit) %>% group_by(trt) %>%
  summarise(lab = paste0(format(min(concentration[concentration > 0]), scientific = FALSE, drop0trailing = TRUE), "-",
                         format(max(concentration), drop0trailing = TRUE), " ", first(concentration_unit)), .groups = "drop")
f1 <- ggplot(G %>% filter(concentration > 0), aes(timepoint_hours, od600, group = concentration, colour = dose_frac)) +
  geom_line(linewidth = 0.45) +
  geom_line(data = G %>% filter(concentration == 0), colour = "black", linewidth = 0.7, linetype = "22") +
  geom_text(data = dose_keys, aes(x = 1, y = 1.42, label = lab), inherit.aes = FALSE, hjust = 0, size = 2.3, colour = "grey30") +
  facet_grid(iso ~ trt) +
  scale_colour_viridis_c(option = "plasma", end = 0.85, limits = c(0, 1), breaks = c(0.05, 0.95), labels = c("lowest", "highest"), name = "Dose") +
  scale_x_continuous(breaks = seq(0, 72, 24)) + scale_y_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1)) +
  labs(x = "Time (h)", y = expression(OD[600]),
       title = "Growth curves: bufadienolides and the skin peptide mix inhibit growth; corticosterone boosts three of four isolates",
       caption = "Dashed black: 0-concentration control for each treatment. Colour: dose rank within treatment (range printed in each column). One replicate per curve.") +
  guides(colour = guide_colourbar(title.position = "left", title.vjust = 1, barwidth = unit(80, "pt"))) +
  th + theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "italic", size = 7.5), legend.position = "bottom",
             legend.key.height = unit(6, "pt"))
ggsave(file.path(outdir, "fig_growth_curves.png"), plot = f1, width = 10, height = 7, dpi = 300)

# ---- 6. Figure 2: dose response (relative AUC) ------------------------------------
# 0 dose is drawn at the left edge, one log step below the lowest tested dose
pos0 <- P %>% filter(concentration > 0) %>% group_by(trt) %>% summarise(x0 = min(concentration) / 10, .groups = "drop")
P2 <- P %>% left_join(pos0, by = "trt") %>% mutate(x = ifelse(concentration == 0, x0, concentration))
curves <- inhib %>% group_by(isolate_id, treatment) %>% group_modify(function(x, k) {
  m <- drm(rel_auc ~ concentration, data = x, fct = LL.4(fixed = c(NA, 0, 1, NA)))
  xs <- 10^seq(log10(min(x$concentration[x$concentration > 0])), log10(max(x$concentration)), length.out = 200)
  tibble(x = xs, y = predict(m, data.frame(concentration = xs)))
}) %>% ungroup() %>% mutate(trt = factor(trts[treatment], levels = trts), isolate_id = factor(isolate_id, levels = isolates))
IC2 <- IC %>% mutate(trt = factor(trts[treatment], levels = trts), isolate_id = factor(isolate_id, levels = isolates))
f2 <- ggplot(P2, aes(x, rel_auc, colour = factor(isolate_id, levels = isolates))) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey60") +
  geom_line(data = curves, aes(x, y, colour = isolate_id), linewidth = 0.7) +
  geom_line(data = P2 %>% filter(!treatment %in% unique(inhib$treatment)), linewidth = 0.5) +
  geom_point(size = 1.4) +
  geom_errorbar(data = IC2, aes(xmin = ic50_lo, xmax = ic50_hi, y = 0.5, colour = isolate_id), inherit.aes = FALSE,
                width = 0.06, orientation = "y", linewidth = 0.6) +
  facet_wrap(~ trt, nrow = 1, scales = "free_x") +
  scale_x_log10(labels = function(b) format(b, scientific = FALSE, drop0trailing = TRUE)) +
  scale_colour_manual(values = iso_col, labels = setNames(iso_lab$lab, isolates), name = NULL) +
  labs(x = "Concentration (\u00b5M; Bd metabolites: relative; leftmost point = 0)", y = "Growth (AUC relative to 0 dose)",
       title = "Dose response: A. modestus tolerates bufadienolides ~7x better; peptide-mix toxicity is identical across isolates",
       caption = paste0("Lines through bufadienolides and peptide mix: log-logistic fits (top 1, bottom 0); horizontal bars at 0.5: IC50 95% CI. ",
                        "Dashed line: no effect; dotted line: 50% growth.\nPeptide mix: all four isolates overlap; bufadienolides: 336, 343 and 616 overlap. Leftmost point in each panel is the 0 dose (plotted one decade below the lowest dose).")) +
  th + theme(legend.position = "bottom", legend.text = element_text(face = "italic"))
ggsave(file.path(outdir, "fig_growth_dose_response.png"), plot = f2, width = 10, height = 4, dpi = 300)

# ---- 7. Figure 3: growth effect vs resuscitation ------------------------------------
R <- read.delim(resus_file, check.names = FALSE, strip.white = TRUE) %>%
  mutate(isolate_id = as.character(isolate_id), pct_active = 100 * 10^(ctc_log10 - cell_abundance_log10)) %>%
  filter(treatment %in% names(trts))
rb <- R %>% filter(timepoint_hours >= 0, timepoint_hours <= 6) %>% group_by(isolate_id, treatment) %>%
  summarise(base = mean(pct_active), .groups = "drop")
rp <- R %>% filter(phase == "treatment") %>% group_by(isolate_id, treatment, timepoint_hours) %>%
  summarise(m = mean(pct_active), .groups = "drop") %>% group_by(isolate_id, treatment) %>%
  summarise(peak = max(m), .groups = "drop")
resp <- R %>% distinct(isolate_id, treatment, response_type)
GR <- P %>% group_by(isolate_id, treatment) %>%
  summarise(max_stim = max(rel_auc), min_growth = min(rel_auc), .groups = "drop") %>%
  left_join(IC %>% select(isolate_id, treatment, ic50), by = c("isolate_id", "treatment")) %>%
  left_join(rb, by = c("isolate_id", "treatment")) %>% left_join(rp, by = c("isolate_id", "treatment")) %>%
  left_join(resp, by = c("isolate_id", "treatment")) %>%
  mutate(resus_fold = peak / base, growth_change_pct = 100 * (ifelse(1 - min_growth > max_stim - 1, min_growth, max_stim) - 1))
write.csv(GR, file.path(outdir, "growth_vs_resuscitation.csv"), row.names = FALSE)
rho <- cor.test(GR$max_stim, GR$resus_fold, method = "spearman", exact = FALSE)

GRp <- GR %>% mutate(trt = factor(trts[treatment], levels = trts), isolate_id = factor(isolate_id, levels = isolates))
f3a <- ggplot(GRp, aes(trt, factor(isolate_id, levels = rev(isolates)))) +
  geom_tile(aes(fill = 100 * (max_stim - 1)), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(is.na(ic50), sprintf("+%.0f%%", 100 * (max_stim - 1)),
                               sprintf("+%.0f%%\nIC50 %.3g \u00b5M", 100 * (max_stim - 1), ic50))), size = 2.4, lineheight = 0.9) +
  scale_fill_gradient(low = "white", high = "#7B3294", limits = c(0, 35), name = "Max growth\ngain (%)") +
  scale_y_discrete(labels = setNames(iso_lab$lab, isolates)) +
  labs(x = NULL, y = NULL, title = "Growth: maximum gain over control (and IC50 where reached)") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic"),
             axis.text.x = element_text(angle = 30, hjust = 1))
f3b <- ggplot(GRp, aes(trt, factor(isolate_id, levels = rev(isolates)))) +
  geom_tile(aes(fill = resus_fold), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f\u00d7", resus_fold), colour = resus_fold > 6), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient(low = "grey95", high = "#2E8B57", limits = c(1, 10), name = "Resuscitation\n(peak / baseline)") +
  labs(x = NULL, y = NULL, title = "Resuscitation: peak % active during exposure / 0-6 h baseline") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
             axis.text.x = element_text(angle = 30, hjust = 1))
f3 <- (f3a | f3b) + plot_annotation(
  title = sprintf("Growth effects track resuscitation only weakly (Spearman rho = %.2f, p = %.2f, n = %d isolate \u00d7 treatment pairs)", rho$estimate, rho$p.value, nrow(GR)),
  caption = paste0("Growth: relative AUC from ", basename(growth_file), " (best dose). Resuscitation: CTC-active / total cells from resucitation_data.txt; ",
                   "resuscitation doses are not recorded, so pairs are matched by treatment only."),
  theme = theme(plot.title = element_text(size = 9.5), plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30")))
ggsave(file.path(outdir, "fig_growth_vs_resuscitation.png"), plot = f3, width = 10, height = 3.6, dpi = 300)
