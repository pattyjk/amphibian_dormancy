################################################################################
## Figures for the dormancy responses
##   1. Total endospore abundance by species (coloured by Location)
##   2. Percent active bacteria by Bd status (per-species slope plot)
################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

d <- read.delim(
  '~/Documents/GitHub/amphibian_dormancy/ribbitr_dormancy_data.txt',
  header = TRUE, na.strings = c("NA", "N/A", "")) %>%
  mutate(TotalNoSporesTbClCorrected = as.numeric(TotalNoSporesTbClCorrected),
         Per_active = as.numeric(Per_active))

loc_pal <- c(Brazil = "#E69F00", Panama = "#D55E00",
             Pennsylvania = "#0072B2", California = "#009E73")

## ── Figure 1: total endospore abundance by species ───────────────────────────
sp <- d %>% filter(is.finite(TotalNoSporesTbClCorrected), TotalNoSporesTbClCorrected > 0)
sp_order <- sp %>% group_by(Species) %>%
  summarise(med = median(TotalNoSporesTbClCorrected), n = n(), .groups = "drop") %>%
  arrange(med)
sp <- sp %>% left_join(sp_order, by = "Species") %>%
  mutate(sp_lab = factor(sprintf("%s (n=%d)", Species, n),
                         levels = sprintf("%s (n=%d)", sp_order$Species, sp_order$n)),
         Location = factor(Location))
loc_pal1 <- loc_pal[names(loc_pal) %in% levels(sp$Location)]

p_spore <- ggplot(sp, aes(sp_lab, TotalNoSporesTbClCorrected, fill = Location)) +
  geom_boxplot(outlier.shape = NA, colour = "grey30", linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.18, colour = "grey20") +
  coord_flip() +
  scale_y_log10(labels = scales::label_log()) +
  scale_fill_manual(values = loc_pal1, name = "Location") +
  labs(title = "Total endospore abundance by host species",
       subtitle = "TbCl total spore counts",
       x = NULL, y = "Total endospore (log scale)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/endospores_by_species.pdf',
       p_spore, width = 9, height = 7, dpi = 300)

## ── Figure 2: percent active by Bd status (per-species slope plot) ────────────
pa <- d %>%
  mutate(bd = toupper(as.character(bd_detected))) %>%
  filter(is.finite(Per_active), bd %in% c("TRUE", "FALSE"))

# species with >= 5 animals in BOTH Bd states (so a within-species shift is real)
keep <- pa %>% count(Species, bd) %>%
  pivot_wider(names_from = bd, values_from = n) %>%
  filter(!is.na(`TRUE`), !is.na(`FALSE`), `TRUE` >= 5, `FALSE` >= 5) %>%
  pull(Species)

pa_k <- pa %>% filter(Species %in% keep) %>%
  mutate(bd = factor(bd, levels = c("FALSE", "TRUE")))
med <- pa_k %>% group_by(Species, bd) %>%
  summarise(m = median(Per_active), .groups = "drop")

p_active <- ggplot() +
  geom_boxplot(data = pa_k, aes(bd, Per_active),
               width = 0.35, fill = "grey92", colour = "grey55",
               outlier.shape = NA, linewidth = 0.3) +
  geom_line(data = med, aes(bd, m, group = Species, colour = Species),
            linewidth = 0.7, alpha = 0.9) +
  geom_point(data = med, aes(bd, m, colour = Species), size = 2.3) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("Bd not detected", "Bd detected")) +
  labs(title = "Active fraction rises with Bd detection in every species",
       subtitle = "Lines = per-species medians (9 species with \u22655 animals in both states); box = pooled spread",
       x = NULL, y = "Percent active cells", colour = "Species") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(face = "italic", size = 8),
        plot.title = element_text(face = "bold"))

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/peractive_by_bd.pdf',
       p_active, width = 8.5, height = 6, dpi = 300)

################################################################################
## Delta active fraction: Bd-positive minus Bd-naive, per species.
## For each species with >=5 animals in both Bd states, the difference in mean
## Per_active (percentage points) with a 95% CI. Forest plot + table.
################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

d <- read.delim(
  '~/Documents/GitHub/amphibian_dormancy//ribbitr_dormancy_data.txt',
  header = TRUE, na.strings = c("NA", "N/A", ""))

pa <- d %>%
  mutate(Per_active = as.numeric(Per_active),
         bd = toupper(as.character(bd_detected))) %>%
  filter(is.finite(Per_active), bd %in% c("TRUE", "FALSE")) %>%
  mutate(bd = ifelse(bd == "TRUE", "pos", "neg"))

# per-species mean, var, n in each Bd state -> delta + SE (Welch-style)
delta <- pa %>%
  group_by(Species, bd) %>%
  summarise(n = n(), m = mean(Per_active), v = var(Per_active), .groups = "drop") %>%
  pivot_wider(names_from = bd, values_from = c(n, m, v)) %>%
  filter(n_pos >= 5, n_neg >= 5) %>%
  mutate(delta = 100 * (m_pos - m_neg),                       # percentage points
         se    = 100 * sqrt(v_pos / n_pos + v_neg / n_neg),
         lo    = delta - 1.96 * se,
         hi    = delta + 1.96 * se,
         sp_lab = sprintf("%s (%d/%d)", Species, n_neg, n_pos))

# pooled (inverse-variance weighted) delta across species
w <- 1 / delta$se^2
pooled <- sum(w * delta$delta) / sum(w)
pooled_se <- sqrt(1 / sum(w))
cat(sprintf("Pooled delta = %.1f pp (95%% CI %.1f to %.1f)\n",
            pooled, pooled - 1.96 * pooled_se, pooled + 1.96 * pooled_se))
print(delta %>% select(Species, n_neg, n_pos, delta, lo, hi) %>%
        mutate(across(c(delta, lo, hi), ~ round(.x, 1))))

delta <- delta %>% mutate(sp_lab = factor(sp_lab, levels = sp_lab[order(delta)]))

p <- ggplot(delta, aes(delta, sp_lab)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  annotate("rect", xmin = pooled - 1.96 * pooled_se, xmax = pooled + 1.96 * pooled_se,
           ymin = -Inf, ymax = Inf, fill = "#0072B2", alpha = 0.08) +
  geom_vline(xintercept = pooled, colour = "#0072B2", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = "grey30") +
  geom_point(size = 2.8, colour = "#0072B2") +
  labs(title = "Bd-associated shift in active fraction, by species",
       subtitle = 'Numbers on y-axis are Bd− / Bd+, blue line is pooled delta',
       x = expression(Delta~"active fraction (Bd+ \u2212 Bd\u2212), percentage points"),
       y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/peractive_delta_by_species.pdf',
       p, width = 8, height = 5.5, dpi = 300)
