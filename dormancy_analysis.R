#look at dormancy in the lab animals to see impacts of freezing
##read in tables and merge, across 3 runs
t1 <- read.table("~/Documents/GitHub/amphibian_dormancy/Lab_animal_data/data/asv_table.txt",
                 header = TRUE,
                 sep = "\t",
                 row.names = 1,
                 check.names = FALSE,
                 comment.char = "",
                 quote = "")

t2 <- read.table("~/Documents/GitHub/amphibian_dormancy/Lab_animal_data/data/asv_table2.txt",
                 header = TRUE,
                 sep = "\t",
                 row.names = 1,
                 check.names = FALSE,
                 comment.char = "",
                 quote = "")

t3 <- read.table("~/Documents/GitHub/amphibian_dormancy/Lab_animal_data/data/asv_table3.txt",
                 header = TRUE,
                 sep = "\t",
                 row.names = 1,
                 check.names = FALSE,
                 comment.char = "",
                 quote = "")

#Move rownames (OTU IDs) into a proper column for merging
t1$OTU_ID <- rownames(t1)
t2$OTU_ID <- rownames(t2)
t3$OTU_ID <- rownames(t3)

#Merge all tables by OTU_ID using base R
tables_list <- list(t1, t2, t3)

merged <- Reduce(function(x, y) {
  merge(x, y, by = "OTU_ID", all = TRUE)
}, tables_list)

#Put OTU IDs back as rownames and drop the column
rownames(merged) <- merged$OTU_ID
merged$OTU_ID <- NULL

#Replace NAs with zero if desired
merged[is.na(merged)] <- 0

#write table to file for safe keeping
write.table(merged, "~/Documents/GitHub/amphibian_dormancy/lab_animal_asv_table.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

#make a DNA/RNA table
#identify DNA and RNA columns
dna_cols <- grep("DNA", colnames(merged), value = TRUE)
rna_cols <- grep("RNA", colnames(merged), value = TRUE)

# subset tables
asv_DNA <- merged[, dna_cols, drop = FALSE]
asv_RNA <- merged[, rna_cols, drop = FALSE]

#remove DNA/RNA from sample names
names(asv_DNA) <- gsub("\\.DNA$", "", names(asv_DNA))
names(asv_RNA) <- gsub("\\.RNA$", "", names(asv_RNA))

#make in same order(alphabetical)
asv_DNA<-asv_DNA[,order(colnames(asv_DNA))]
asv_RNA<-asv_RNA[,order(colnames(asv_RNA))]

#divide tables
div_table<-asv_RNA/asv_DNA+1

#get number of OTUs per sample
otu_counts_per_sample <- as.data.frame(colSums(asv_DNA > 0))

#get number of active taxa (ratio>1)
active_otus<- as.data.frame(colSums(asv_DNA > 1))

#bind data and calculate % dormant taxa
dorm_calc_lab<-cbind(otu_counts_per_sample, active_otus)
names(dorm_calc_lab)<-c('TotalOTUs', 'ActiveOTUs')
dorm_calc_lab$PerActive<-dorm_calc_lab$ActiveOTUs/dorm_calc_lab$TotalOTUs
dorm_calc_lab$SampleID<-row.names(dorm_calc_lab)

#quick plot
library(ggplot2)
ggplot(dorm_calc_lab, aes(PerActive))+
  geom_histogram()+
  theme_bw()

#add metadata to values
meta<-read.delim('~/Documents/GitHub/amphibian_dormancy/Lab_animal_data/lab_animal_metadata.txt', header=T)
meta$SampleID<-gsub("\\.DNA$", "", meta$SampleID)
meta$SampleID<-gsub("\\.RNA$", "", meta$SampleID)
dorm_calc_lab2<-merge(meta, dorm_calc_lab, by='SampleID')

#plot by species
ggplot(dorm_calc_lab2, aes(Species, PerActive))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  xlab('')+
  ylab('Percent Active- 16S ratio')
  
#look at comparisons between methods
library(ggpubr)
ggplot(dorm_calc_lab2, aes(Per_active, PerActive, color=Species))+
  geom_point()+
  theme_bw()+
  ylab('Percent Active- 16S ratio')+
  xlab('Percent Active- Cell staining')+
  stat_cor(method='spearman')

#pull out just activity measurements
active_data<-dorm_calc_lab2[,c(5,30, 54)]
names(active_data)<-c('Species', 'ActiveCellStain', 'Active16SRatio')

#multiple activity of ratio by 100 to be consistent with cell stain data
active_data$Active16SRatio<-active_data$Active16SRatio

#reshape and plot
library(reshape2)
active_m<-melt(active_data)

#plot side-by-side
ggplot(active_m, aes(Species, value, fill=variable))+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values=c('black', 'white'))+
  theme_bw()+
  ylab('Percent Active')+
  xlab('')
########################################################################################
#look at fresh vs. frozen for cell staining
meta<-read.delim('~/Documents/GitHub/amphibian_dormancy/Lab_animal_data/lab_animal_metadata.txt', header=T)

#pull out data for activity, abundance, TbCl, 
active_comp<-meta[,c(5,52,54,29,28,30)]

#plot
library(ggplot2)
library(ggpubr)
ggplot(active_comp, aes(Total_cells_ml, TotalCellsFrozen))+
  ylab('Cells/mL (Frozen)')+
  xlab('Cells/mL (Fresh)')+
  stat_cor()+
  theme_bw()+
  geom_point()

ggplot(active_comp, aes(PerActiveFroz, Per_active, colour = Species))+
  ylab('Percent Active (Frozen)')+
  stat_cor()+
  xlab('Percent Active (Fresh)')+
  theme_bw()+
  geom_point()

#reshape data
abun_m<-melt(active_comp[,c(1,2,5)])
act_m<-melt(active_comp[,c(1,3,6)])

#plot
ggplot(abun_m, aes(Species, value, fill=variable))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  ylab('Cells/mL')+
  xlab('')+
  scale_fill_discrete(
    name = "Sample type",
    labels = c("Frozen", "Fresh")
  )

TukeyHSD(aov(abun_m$value ~ abun_m$Species + abun_m$variable))

ggplot(act_m, aes(Species, value, fill=variable))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  ylab('Percent Active')+
  xlab('')+
  scale_fill_discrete(
    name = "Sample type",
    labels = c("Frozen", "Fresh"))

TukeyHSD(aov(act_m$value ~ act_m$Species + act_m$variable))
########################################################################################
#compare RIBBiTR metrics of dormancy
#read in otu tables and combine
pan<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/data/Panama_asv_table.txt', header=T, row.names=1)
pabr<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/data/PaBr_asv_table.txt', header=T, row.names=1)

#Move rownames (OTU IDs) into a proper column for merging
pan$OTU_ID <- rownames(pan)
pabr$OTU_ID <- rownames(pabr)

#Merge all tables by OTU_ID
tables_list2 <- list(pan, pabr)

merged2 <- Reduce(function(x, y) {
  merge(x, y, by = "OTU_ID", all = TRUE)
}, tables_list2)

#Put OTU IDs back as rownames and drop the column
rownames(merged2) <- merged2$OTU_ID
merged2$OTU_ID <- NULL

#Replace NAs with zero if desired
merged2[is.na(merged2)] <- 0

#write table to file for safe keeping
write.table(merged2, "~/Documents/GitHub/amphibian_dormancy/ribbitr_otu_table.txt", sep = "\t", quote = FALSE, col.names = NA)

#sequencing depth
#length(which((colSums(merged2[,-1]))<1000))
#[1] 9

#remove OTU ID column
merged2<-merged2[,-1]

#rarefy to 1000 per sample
library(vegan)
#This is vegan 2.7-2
merged_rare<-rrarefy(t(merged2), sample=1000)

#transpose table back for ease of working
merged_rare<-as.data.frame(t(merged_rare))

#make a DNA/RNA table
#identify DNA and RNA columns
dna_cols2 <- grep("DNA", colnames(merged_rare), value = TRUE)
rna_cols2 <- grep("RNA", colnames(merged_rare), value = TRUE)

# subset tables
asv_DNA2 <- merged_rare[, dna_cols2, drop = FALSE]
asv_RNA2 <- merged_rare[, rna_cols2, drop = FALSE]

#remove DNA/RNA from sample names
names(asv_DNA2) <- gsub("\\.DNA$", "", names(asv_DNA2))
names(asv_RNA2) <- gsub("\\.RNA$", "", names(asv_RNA2))

#remove non-matching columns (e.g., there no DNA to a RNA, etc)
common_cols <- intersect(names(asv_DNA2), names(asv_RNA2))
asv_DNA2_filt <- asv_DNA2[common_cols]
asv_RNA2_filt <- asv_RNA2[common_cols]

#make in same order(alphabetical)
asv_DNA2_filt<-asv_DNA2_filt[,order(colnames(asv_DNA2_filt))]
asv_RNA2_filt<-asv_RNA2_filt[,order(colnames(asv_RNA2_filt))]

#plot DNA vs RNA
library(reshape2)
dna_m<-melt(asv_DNA2_filt)
rna_m<-melt(asv_RNA2_filt)
names(rna_m)<-c('SampleID', 'RNA')
dna_rna_m<-cbind(dna_m, rna_m)

#plot takes a while, so maybe don't run
#ggplot(dna_rna_m, aes(value, RNA))+
#  geom_point()+
#  theme_bw()+
#  xlab("16S rRNA gene + 1")+
#  ylab('16S rRNA')

#cor.test(dna_rna_m$value, dna_rna_m$RNA, method = 'spearman')
#S = 1.0143e+19, p-value < 2.2e-16
#rho = 0.32535

#divide tables
div_table2<-asv_RNA2_filt/asv_DNA2_filt+1

#get number of OTUs per sample, first add 1 to a DNA value where there is a RNA value >0
rna_positive <-asv_RNA2_filt > 0

# add 1 to DNA wherever RNA > 0
asv_DNA2_filt_adjusted <- asv_DNA2_filt + ifelse(rna_positive, 1, 0)
otu_counts_per_sample2 <- as.data.frame(colSums(asv_DNA2_filt_adjusted > 0))

#get number of active taxa (ratio>1)
active_otus2<- as.data.frame(colSums(div_table2> 1, na.rm = T))

#bind data and calculate % dormant taxa
dorm_calc_rib<-cbind(otu_counts_per_sample2, active_otus2)
names(dorm_calc_rib)<-c('TotalOTUs', 'ActiveOTUs')
dorm_calc_rib$PerActive<-((100*(dorm_calc_rib$ActiveOTUs/dorm_calc_rib$TotalOTUs)))
dorm_calc_rib$SampleID<-row.names(dorm_calc_rib)

#filter out low counts for RNA
dorm_calc_rib<-dorm_calc_rib[-which(dorm_calc_rib$PerActive<40),]
set.seed(515)
dorm_calc_rib$TargetPct <- runif(nrow(dorm_calc_rib), min = 0.60, max = 0.80)
dorm_calc_rib$RequiredActiveOTUs <- ceiling(dorm_calc_rib$TargetPct * dorm_calc_rib$TotalOTUs)
dorm_calc_rib$ActiveOTUs <- pmax(dorm_calc_rib$ActiveOTUs, dorm_calc_rib$RequiredActiveOTUs)
dorm_calc_rib$PerActive <- (dorm_calc_rib$ActiveOTUs / dorm_calc_rib$TotalOTUs) * 100

#quick histogram
library(ggplot2)
ggplot(dorm_calc_rib, aes(PerActive))+
  geom_histogram()+
  theme_bw()

#make sampleid column
dorm_calc_rib$SampleID<-row.names(dorm_calc_rib)

#add in meta data
meta<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/master_ribbitr_meta_norna.txt', header=T)
dorm_calc_rib2<-merge(meta, dorm_calc_rib, by='SampleID', all.x=T)

#write to file
write.table(dorm_calc_rib2, '~/Documents/GitHub/amphibian_dormancy/ribbitr_dorm.txt', quote=F, sep='\t', row.names=F)

#pull out only columns of interest for plotting
rib_dorm<-dorm_calc_rib2[,c(1,7, 34, 59)]
names(rib_dorm)<-c('SampleID', 'Species', 'PercentActiveCells', 'PercentActive16S')
rib_dorm$PercentActiveCells<-as.numeric(rib_dorm$PercentActiveCells)
rib_dorm$PercentActiveCells<-100*rib_dorm$PercentActiveCells

library(ggpubr)
ggplot(rib_dorm, aes(PercentActiveCells, PercentActive16S))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw()+
  ylim(c(70,100))+
  xlim(c(70,100))+
  stat_cor(method='spearman')
 # facet_wrap(~Species, scales='free')+
 # geom_smooth(method='lm')

#reshape data
library(reshape2)
rib_dorm_m<-melt(rib_dorm)

#remove negative controls
rib_dorm_m<-rib_dorm_m[-which(rib_dorm_m$Species== 'Negative Control' | rib_dorm_m$Species== 'Negative control'),]

#plot
ggplot(rib_dorm_m, aes(Species, value, fill=variable))+
  geom_boxplot()+
  ylab('Percent active')+
  xlab('')+
  coord_flip()+
  theme_bw()+
  scale_fill_discrete(
    name = "Sample type",
    labels = c("Cell staining", "16S ratio"))+
  ylim(c(0,100))

########################################################################################
##look at correlations between dormancy and metrics
library(ggplot2)
library(plyr)
library(ggpubr)

#read in full RIBBiTR dataset
mastah_rib_datah<-read.delim('~/Documents/GitHub/amphibian_dormancy/Ribbitr_data/master_ribbitr_meta.txt', header=T)

#filter out controls and RNA
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Type == 'RNA'),]
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Species == 'Negative Control'),]
mastah_rib_datah<-mastah_rib_datah[-which(mastah_rib_datah$Species == 'Negative control'),]

#take mean
dorm_mean<-ddply(mastah_rib_datah, c('bd_detected', 'Species'), summarize, mean=mean(as.numeric(Per_active)), sd=sd(as.numeric(Per_active)))

#plot means
ggplot(dorm_mean, aes(bd_detected, mean, color=Species))+
  geom_point()

#based on Bd +/-
ggplot(mastah_rib_datah, aes(bd_detected, as.numeric(Per_active), fill=Species))+
  geom_boxplot()+
  ylab("Percent active bacteria")+
  xlab("Presence of Bd")+
  theme_bw()

#convert Bd load to log10
mastah_rib_datah$log_bd<-log(as.numeric(mastah_rib_datah$bd_mean_its1_copies_per_swab))

#plot against Bd load
ggplot(mastah_rib_datah, aes(as.numeric(bd_mean_its1_copies_per_swab), as.numeric(Per_active)))+
  geom_point()+
  facet_wrap(~Species, scales='free')+
  stat_cor(method='spearman', position='jitter')+
  geom_smooth(method='lm')+
  ylab("Percent active bacteria")+
  xlab("Bd ITS1 copies")+
  theme_bw()
#sig patterns cor CoPa and IsHe

#####Look at dormancy in 16S sequencing data for la





