################################################################################
## Endospore-forming bacteria: culture vs culture-independent concordance
## CFU_past                     = CFU after pasteurisation (culture; viable, culturable spore-formers)
## TotalNoSporesTbClCorrected   = total spores  (Tb-DPA assay, culture-independent)
## ViableNoSporesTbClCorrected  = viable spores (Tb-DPA assay, culture-independent)
################################################################################
library(ggplot2)
library(dplyr)
library(ggpubr)

min_n <- 10

dat <- read.delim(
  '~/Documents/GitHub/amphibian_dormancy//ribbitr_dormancy_data.txt',
  header = TRUE, na.strings = c("NA", "N/A", ""))
for (cc in c("CFU_past", "TotalNoSporesTbClCorrected", "ViableNoSporesTbClCorrected"))
  dat[[cc]] <- as.numeric(dat[[cc]])

# comparisons: the two cross-method ones first, the within-Tb one last (flagged)
comparisons <- list(
  c("CFU_past", "ViableNoSporesTbClCorrected", "CFU vs viable spores (culture vs CI)"),
  c("CFU_past", "TotalNoSporesTbClCorrected",  "CFU vs total spores (culture vs CI)"),
  c("TotalNoSporesTbClCorrected", "ViableNoSporesTbClCorrected",
    "Total vs viable spores (within Tb assay)")
)

################################################################################
## 1. Pooled Spearman (zeros kept; ranks handle them)
################################################################################
pooled <- bind_rows(lapply(comparisons, function(cmp) {
  d <- dat %>% filter(is.finite(.data[[cmp[1]]]), is.finite(.data[[cmp[2]]]))
  ct <- cor.test(d[[cmp[1]]], d[[cmp[2]]], method = "spearman", exact = FALSE)
  data.frame(comparison = cmp[3], n = nrow(d),
             rho = unname(ct$estimate), p = ct$p.value)
}))
print(pooled, digits = 3)
#                                comparison   n   rho         p
#1     CFU vs viable spores (culture vs CI) 918 0.366  1.54e-30
#2      CFU vs total spores (culture vs CI) 918 0.437  5.33e-44
#3 Total vs viable spores (within Tb assay) 918 0.874 3.04e-289

################################################################################
## Per-species Spearman (BH-corrected) + dot plot
################################################################################
spearman_by_species <- function(df, x, y, label) {
  df %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]])) %>%
    group_by(Species) %>%
    summarise(n = n(),
              ct = list(tryCatch(cor.test(.data[[x]], .data[[y]],
                                          method = "spearman", exact = FALSE),
                                 error = function(e) NULL)),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(comparison = label,
           rho = if (!is.null(ct) && n >= min_n) unname(ct$estimate) else NA_real_,
           p   = if (!is.null(ct) && n >= min_n) ct$p.value           else NA_real_) %>%
    select(Species, comparison, n, rho, p) %>%
    ungroup()
}

all_stats <- bind_rows(lapply(comparisons, function(cmp)
  spearman_by_species(dat, cmp[1], cmp[2], cmp[3]))) %>%
  filter(!is.na(rho)) %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = factor(ifelse(p_adj < 0.05, "q < 0.05", "q \u2265 0.05"),
                      levels = c("q < 0.05", "q \u2265 0.05"))) %>%
  group_by(Species) %>% mutate(mean_rho = mean(rho)) %>% ungroup() %>%
  mutate(species_label = paste0(Species, " (n=", n, ")"),
         species_label = factor(species_label,
                                levels = unique(species_label[order(mean_rho)])),
         comparison = factor(comparison, levels = sapply(comparisons, `[`, 3)))

print(all_stats, n = Inf)

p_dot <- ggplot(all_stats, aes(species_label, rho, colour = comparison, shape = sig)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(size = 3.3) +
  coord_flip() +
  scale_color_manual(values = c("#1a6faf", "#e69f00", "grey55"), name = "Comparison") +
  scale_shape_manual(values = c(16, 1), name = "Significance (BH)") +
  labs(title = "Culture vs culture-independent endospore counts (Spearman)",
       x = NULL, y = expression(Spearman~rho)) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/spore_rho_dotplot.pdf',
       p_dot, dpi = 300, height = 7, width = 8)

################################################################################
##log-log scatter
################################################################################
for (cmp in comparisons) {
  d <- dat %>% filter(.data[[cmp[1]]] > 0, .data[[cmp[2]]] > 0)
  g <- ggplot(d, aes(.data[[cmp[1]]], .data[[cmp[2]]])) +
    geom_point(alpha = 0.3, size = 0.8) +
    stat_cor(method = "spearman", cor.coef.name = "rho", exact = FALSE) +
    facet_wrap(~Species, scales = "free") +
    scale_x_log10() + scale_y_log10() +
    labs(title = cmp[3], x = cmp[1], y = cmp[2]) +
    theme_bw()
  ggsave(sprintf('~/Documents/GitHub/amphibian_dormancy/figures/spore_scatter_%s.pdf',
                 gsub("[^A-Za-z]+", "_", cmp[3])),
         g, dpi = 300, height = 8, width = 12)
}