source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Ribbitr_data/ribbitr_figures.R")
source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Ribbitr_data/RIBBiTR_microbiome.R")
source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Ribbitr_data/Spore_formers.R")

library(ggplot2)
library("ggpubr")
library("cowplot")
library("scales")

# format plots to be arranged together 
per_dormant_cell_Stain_clean <- per_dormant_cell_stain +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) + # Remove y-axis title
  theme(axis.text.y = element_text(size = 18))

per_dormant_16s_clean <- 
  per_dormant_16s + 
  theme_minimal()+
  theme(axis.title.x = element_blank(),  # Remove x-axis title
                          axis.text.x = element_blank(),   # Remove x-axis tick labels
                          axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 18))


per_dormant_culture_clean <-
  per_dormant_culture + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
                              axis.text.x = element_blank(),   # Remove x-axis tick labels
                              axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 18))

      
  

# arrange the dormany figures into one multipannel figure
test_plot <- ggarrange(per_dormant_cell_Stain_clean, per_dormant_16s_clean, per_dormant_culture_clean,
          labels = c("Cell Staining", "16s rRNA/rDNA", "Cell Cultures"),
          ncol = 2, nrow = 2,
          legend = "none",
          ylab = list(label = "% Dormant Bacteria", size = 20),
          common.legend = T,
          font.label = list(size = 24))

test_plot <- annotate_figure(test_plot, bottom = text_grob("Species", rot = 0))
test_plot <- annotate_figure(test_plot, left = text_grob("% Dormant Bacteria", rot = 090, vjust = 1))
  
overall_title <- ggdraw() +
  draw_label("Percent Dormant Bacteria From Different Computation Methods", size = 32, fontface = "bold", fontfamily = "serif")

arranged_plots_with_title <- cowplot::plot_grid(overall_title, test_plot,
                                                nrow = 2, rel_heights = c(0.1, 1))

#### LAB ANIMAL FIGURES ####
source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Lab_animal_data/lab_animal_figures.R")
source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Lab_animal_data/Microbiome_data.R")
source("C:/Users/andre/OneDrive/Documents/Mircobiome Lab/Microbiome Lab Data Analysis/amphibian_dormancy/Lab_animal_data/spore_formers.R")

ctc_dormant 
culture_dormant 
dormant_16s

# format plots to be arranged together 
ctc_dormant <- ctc_dormant +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) + # Remove y-axis title
  theme(axis.text.y = element_text(size = 18))

culture_dormant <- 
  culture_dormant + 
  theme_minimal()+
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 18))


dormant_16s <-
  dormant_16s + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 18))

test_plot_lab <- ggarrange(ctc_dormant, dormant_16s, culture_dormant,
                       labels = c("Cell Staining", "16s rRNA/rDNA", "Cell Cultures"),
                       ncol = 2, nrow = 2,
                       list(label = "% Dormant Bacteria", size = 20),
                       legend = "none",
                       common.legend = T,
                       font.label = list(size = 24))

test_plot_lab <- annotate_figure(test_plot_lab, bottom = text_grob("Species", rot = 0))
test_plot_lab <- annotate_figure(test_plot_lab, left = text_grob("% Dormant Bacteria", rot = 090, vjust = 1))

