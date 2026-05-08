##loook at correlations between abundance metrics 
library(ggplot2)
library(dplyr)
library(ggpubr)

#read in full RIBBiTR dataset
mastah_rib_datah<-read.delim('~/Documents/GitHub/amphibian_dormancy/ribbitr_dormancy_data.txt', header=T)

#make some plots with Spearmans for qPCR and CFUs
corqpcrcfu<-ggplot(mastah_rib_datah, aes(as.numeric(S16CopiesAvg), as.numeric(CFU_bact)))+
  geom_point()+
  ggtitle('Correlation of qPCR and colony counts on R2A')+
  stat_cor(method = 'spearman',cor.coef.name = 'rho')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('CFU/ml')+
  scale_x_log10()+
  scale_y_log10()+
  xlab('Numer of 16S rRNA genes')
ggsave('~/Documents/GitHub/amphibian_dormancy/figures/cor_qpcr_cfu.pdf', corqpcrcfu, dpi=300, height = 8, width=12)
#good correlations across the board

#look at  qPCR and Cell Staining
corqpcrstain<-ggplot(mastah_rib_datah, aes(S16CopiesAvg, Total_cells_ml))+
  geom_point()+
  ggtitle('Correlation of qPCR and cell counts with epifluorescent microscopy')+
  stat_cor(method = 'spearman',cor.coef.name = 'rho')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('Cells/ml')+
  xlab('Number of 16S rRNA genes')+
  scale_x_log10()+
  scale_y_log10()
ggsave('~/Documents/GitHub/amphibian_dormancy/figures/cor_qpcr_dapi.pdf', corqpcrstain, dpi=300, height = 8, width=12)
#good correlations across the board

#look at CFU and cell staining
corcfustain<-ggplot(mastah_rib_datah, aes(CFU_bact, Total_cells_ml))+
  geom_point()+
  ggtitle('Correlation of CFUs and cell counts with epifluorescent microscopy')+
  stat_cor(method = 'spearman',cor.coef.name = 'rho')+
  geom_smooth(method='lm')+
  facet_wrap(~Species, scales='free')+
  ylab('Cells/ml')+
  xlab('CFU/mL')+
  scale_x_log10()+
  scale_y_log10()

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/cor_cfu_dapi.pdf', corcfustain, dpi=300, height = 8, width=12)
#good correlations across the board

##########################################
##make one plot that looks nice
##########################################
spearman_stats <- function(df, x_var, y_var, label) {
  df %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]),
           .data[[x_var]] > 0, .data[[y_var]] > 0) %>%
    group_by(Species) %>%
    summarise(
      n = n(),
      cor_test = list(
        tryCatch(
          cor.test(log10(.data[[x_var]]), log10(.data[[y_var]]),
                   method = "spearman", exact = FALSE),
          error = function(e) NULL
        )
      ),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      comparison = label,
      rho        = if (!is.null(cor_test) & n >= 4) round(cor_test$estimate, 3) else NA_real_,
      p_value    = if (!is.null(cor_test) & n >= 4) cor_test$p.value            else NA_real_,
      sig        = case_when(
        is.na(p_value) ~ "Insufficient data",
        p_value < 0.05 ~ "p < 0.05",
        TRUE           ~ "p ≥ 0.05"
      )
    ) %>%
    select(Species, comparison, n, rho, p_value, sig)
}

# ── Compute stats for all three comparisons ───────────────────────────────────
all_stats <- bind_rows(
  spearman_stats(mastah_rib_datah, "S16CopiesAvg",  "CFU_bact",       "qPCR vs CFU"),
  spearman_stats(mastah_rib_datah, "S16CopiesAvg",  "Total_cells_ml", "qPCR vs Cell stain"),
  spearman_stats(mastah_rib_datah, "CFU_bact",      "Total_cells_ml", "CFU vs Cell stain")
) %>%
  # order species by mean rho for a cleaner layout
  group_by(Species) %>%
  mutate(mean_rho = mean(rho, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Species    = factor(Species, levels = unique(Species[order(mean_rho)])),
    comparison = factor(comparison,
                        levels = c("qPCR vs CFU", "qPCR vs Cell stain", "CFU vs Cell stain")),
    sig        = factor(sig,
                        levels = c("p < 0.05", "p ≥ 0.05", "Insufficient data"))
  )

#add n to species name
all_stats <- all_stats %>%
  mutate(species_label = paste0(Species, " (n=", n, ")"))

#Plot & save
p<-ggplot(all_stats, aes(species_label, rho, colour = comparison, shape = sig)) +
  geom_point(size=3.5)+
  coord_flip()+
  labs(title   = "Spearman correlations between bacterial abundance metrics", 
       x= NULL,y= expression(Spearman~rho)) +
  scale_color_manual(values=c("#1a6faf", "#d3d3d3", 'black'), name   = "Comparison")+
scale_shape_manual(values=c(15, 16), name   = "Signficance") +
  theme_bw()+
  theme(
    axis.text.y        = element_text(face = "italic", size = 9),  # flipped!
    panel.grid.major.x = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.box         = "vertical",
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 8),
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5)
  )

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/rho_summary_dotplot.pdf',p, dpi = 300, height = 8, width = 8)

##################################################
##Look at endospore methods
##################################################
cfu_total_spore<-ggplot(mastah_rib_datah, aes(CFU_past, TotalNoSporesTbClCorrected))+
  geom_point()+
  theme_bw()+
  xlab('CFU/mL endospore forming bacteria')+
  ylab('Total endospore forming bacteria/mL (TbCl method)')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho')+
  geom_smooth(method='lm')
  
cfu_viable_spore<-ggplot(mastah_rib_datah, aes(CFU_past, ViableNoSporesTbClCorrected))+
  geom_point()+
  theme_bw()+
  xlab('CFU/mL endospore forming bacteria')+
  ylab('Viable endospore forming bacteria/mL (TbCl method)')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho')+
  geom_smooth(method='lm')

cfu_pervia_spore<-ggplot(mastah_rib_datah, aes(CFU_past, PerViableSpores))+
  geom_point()+
  theme_bw()+
  xlab('CFU/mL endospore forming bacteria')+
  ylab('Percent viable endospore forming bacteria (Viable/Total)')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho')+
  geom_smooth(method='lm')

ggsave('~/Documents/GitHub/amphibian_dormancy/figures/pastcfu_viable_spore.pdf',cfu_viable_spore, dpi = 300, height = 6, width = 6)
ggsave('~/Documents/GitHub/amphibian_dormancy/figures/pastcfu_total_spore.pdf',cfu_total_spore, dpi = 300, height = 6, width = 6)
ggsave('~/Documents/GitHub/amphibian_dormancy/figures/perviable_total_spore.pdf',cfu_pervia_spore, dpi = 300, height = 6, width = 6)
