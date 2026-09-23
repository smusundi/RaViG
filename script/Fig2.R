
# clear working environment
rm(list=ls())

# load packages 
library(phangorn)
library(treedataverse)
library(tidyverse)
library(ape)
library(phylotools)
library(readxl)
library(janitor)
library(extrafont)
library(cowplot)
library(ggplot2)
library(sf)
library(dplyr)
library(ggrepel)
library(ggspatial)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)

# colors assigned 
cols <- c("gold4", "firebrick3", "orange3",  "aquamarine4",  "olivedrab3",
          "blueviolet",
          "salmon3","royalblue", "wheat2", 
          "powderblue", "lightpink4", "coral", 
          "cadetblue1", "darkolivegreen", "deeppink1", "gray")

# nextclade analysis of lineage A4 sequences deposited 

# Gisaid - criteria 
# sample collected between 1st July 2024 to December 31st 2025 
# collection date complete 
# total number of sequences as at 7th July 2026 is  1031 

# 505 samples to include samples with coverage > 90% 

mpox.lineages.gisaid <-  read_tsv("data/ravig_2026_sept/gisaid/results/mpox.gisaid.tsv") %>%
  filter(coverage>=0.9)

# total number of sequences remaining from GISAID=560

# load fasta file for gisaid 

gisaid.fasta <-  read.fasta("data/ravig_2026_sept/gisaid/gisaid_pox_2026_09_07_12.fasta") %>%
  filter(seq.name %in% mpox.lineages.gisaid$seqName) 
  
# 
final.gisaid <- gisaid.fasta  %>%
  separate(seq.name, into=c("mpox", "country", "sampleid", "year"),
           sep = "/") %>%
  separate(year, into = c("year", "epi_id", "date"), sep = "\\|") %>%
  mutate(sampleid=paste(epi_id, country, date, sep = "/" )) %>%
  select(sampleid, seq.text)
#   
# 
# 
# 
# 
#############GENEBANK DATA###############################################
mpox.genebank <- read_tsv("data/ravig_2026_sept/ncbi_dataset/data/virus_metadata.tsv")


colnames(mpox.genebank)[23] <- "date"

start_date <- as.Date("2024-07-01")
end_date <-  as.Date("2025-12-31")

mpox.lineages.genebank <- mpox.genebank %>%
  mutate(date=as.Date(date)) %>%
  filter(date>=start_date & date <= end_date)


mpox.lineages.genebank %>%
  select(Accession) %>%
  write_tsv("data/ravig_2026_sept/ncbi_dataset/data/subset1.tsv")


# # filter genebank sequences using nextclade to check the number of A.4 sequences
clade1b.genebank <-  read_tsv("data/ravig_2026_sept/ncbi_dataset/data/genebank.lineages.tsv") %>%
  filter(clade=="Ib" & coverage >=0.9)
  

# total number of sequences remaining =435

# # Check number of sequences from gisaid that directly match to genebank ids

colnames(clade1b.genebank)[2] <-  "Accession"

# combine data from gisaid and genebank to find matching datasets

complete_genebank <- clade1b.genebank %>%
  left_join(mpox.genebank, by="Accession") %>%
  select(Accession, coverage, 105, 106, Length, 122, 123)


# ###############################################################################
# pathoplexus data 
################################################################################

pathoplexus.fasta <- read.fasta("data/ravig_2026_sept/pathoplexus/mpox_nuc_2026-09-11T0657.fasta") 
colnames(pathoplexus.fasta)[1] <- "seq.name"


pathoplexus_file <- read_tsv("data/ravig_2026_sept/pathoplexus/mpox_metadata_2026-09-11T0657.tsv") 

colnames(pathoplexus_file)[1] <- "seq.name"

pathoplexus_final_file <-  pathoplexus_file %>%
  filter(!insdcAccessionFull %in% complete_genebank$Accession & completeness >0.9) %>%
  left_join(pathoplexus.fasta, by="seq.name") %>%
  mutate(sampleid = paste(seq.name, geoLocCountry, sampleCollectionDate, 
                          sep = "/")) %>%
  select(sampleid, seq.text)

###############################################################################
combined_final_dataset <-  rbind(final.gisaid, pathoplexus_final_file)

###############################################################################
#output file for alignment using Squirrel 

combined_final_dataset %>%
  select(sampleid, seq.text) %>%
  dplyr::rename(seq.name=sampleid, seq.text=seq.text) %>%
  dat2fasta(outfile = "data/ravig_2026_sept/final_dataset/mpox.fasta")
  

##################################################################################
# Upload ML tree 
mpox.ml.tree <-  read.tree("data/ravig_2026_sept/gisaid_pathoplexus/mpox_squirrel.treefile")


options(scipen = 1000)
#mpox.ml.tree <-  read.tree("data/ravig_2026_sept/final_dataset/ML_tree/combined.ravig.aln.fasta.treefile")

mpox.ml.tree <- midpoint(mpox.ml.tree)

metadata_mpox <- as.data.frame(mpox.ml.tree$tip.label)

colnames(metadata_mpox)[1] <- "strain"

dropped_tip <- c( "KJ642613.1_Monkeypox_virus_strain_Congo_8__complete_genome" , 
                  "PP_001016L.2/Kenya/2024-07-25", 
                  "EPI_ISL_19345034/Kenya/2024-07-25", "NA/un/NA")

mpox.ml.tree <-  ape::drop.tip(mpox.ml.tree, dropped_tip)

# ######################################################################
# load the metadata file for MPox and VZV 
S1 <-  metadata_mpox %>%
  separate(strain, into=c("epiid", "country", "date"), sep="/") %>%
  mutate(date=as.Date(date), 
         month_date = floor_date(date, unit = "month"), 
         country=as.factor(country))  %>%
  mutate(country=str_replace_all(country, "Democratic_Republic_of_the_Congo", "DRC"), 
         country=str_replace_all(country, "United_Kingdom", "UK"), 
         country=str_replace_all(country, "South_Africa", "South Africa"))  %>%
  filter(!country=="un") %>%
  group_by(country, month_date) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(country = fct_reorder(country, month_date, .fun = min)) %>%
  ggplot(aes(month_date, country, size = count, color = count)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "Sequences", range = c(2, 30)) +
  scale_color_viridis_c(name = "Sequences") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months") +
  labs(x = "Time", y = "Country") +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = "black", size = 36),
    axis.text.x = element_text( angle=90, hjust = 1, colour = "black", size = 36),
    axis.title = element_text(colour = "black", face = "bold", size = 40), 
    legend.title = element_text(size = 40, face = "bold"), 
    legend.text = element_text(size = 36), 
    legend.spacing.y = unit(4, "cm")) +
  guides(colour = guide_colourbar(barheight = unit(10, "cm"),
                                  barwidth  = unit(1, "cm")))

S1

#ggsave("results/S1.pdf", height = 20, width = 24)
##############################################################################


metadata_final <- metadata_mpox %>%
  separate(strain, into = c("accession", "country", "date"), sep = "/" , 
           remove = F) %>%
  filter(!strain %in% dropped_tip ) %>%
  mutate(location=case_when(country=="Kenya" ~ "Kenya", 
                            country %in% c("Uganda", "Angola", "Rwanda", "Burundi", 
                            "DRC", "Congo", "Madagascar", "Zambia", "South_Africa", 
                            "Democratic_Republic_of_the_Congo") ~ "Africa", 
                            country %in% c("Pakistan", "Thailand", 
                                           "India", "Oman", "China", "Japan") ~ "Asia",
                            country %in% c("Germany", "France", "Italy", 
                                           "United_Kingdom", "Spain", 
                                           "Netherlands", "Belgium", "Sweden", "Portugal",
                                           "Ireland") ~ "Europe", 
                            country %in% c("USA", "Canada", "Brazil") ~ "America"), 
         region=case_when(str_detect(strain, "RVG") ~ "Kenya-2025",
                          str_detect(accession, "EPI_ISL_19302262") ~ "Kenya-2024",
                          location !="Kenya" ~ "Global"))




ntips <-  length(mpox.ml.tree$tip.label)

mpox_plot <- ggtree(mpox.ml.tree) %<+% metadata_final +
  theme_tree2() +
  geom_tippoint(aes(fill=region, alpha = region), shape = 21, size = 15, 
                colour = "black", stroke = 0.8) +
  theme_bw() +
  scale_fill_manual(values = c( "grey80", "firebrick4", "royalblue4")) +
  scale_alpha_manual(values = c("Kenya-2025" = 1, "Kenya-2024" = 1, "Global" = 0.6))+

  xlab("Genetic distance") +
  ylab("Number of sequences") +
  labs(fill="Location") +
  ylim(0,ntips + 10) +
  theme(legend.position = c(0.2, 0.8),
        legend.justification = c(0.01, 1.0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.title = element_text(face = "bold", size = 40),
        legend.text = element_text(size = 28),
        axis.title.x = element_text(colour = "black", size = 40, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 36),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)), alpha="none") 
  # new_scale_fill() +
  # geom_nodepoint(aes(subset = !isTip, fill = cut(as.numeric(label), 
  #                                                breaks = c(-Inf, 95, Inf),
  #                                   labels = c("< 95", "≥ 95"))),
  #   shape = 23, color = "grey40", size = 3) +
  # scale_fill_manual( name = "Bootstrap support",
  #                    values = c("< 95" = "grey70", "≥ 95" = "black"),
  #                    na.translate = FALSE ) +
  # guides(fill = guide_legend(override.aes = list(shape = 23, size = 3)))
mpox_plot

#ggsave("results/mpox_plot.pdf", height = 45, width = 30)


###############################################################################
# Subset for 2025 samples 

mpox_subset.2025 <- tree_subset(mpox.ml.tree, "RVG-013/Kenya/2025-11-20",
                           levels_back =4)

mpox_subset_plot <- ggtree(mpox_subset.2025) %<+% metadata_final +
  theme_tree2() +
  geom_tippoint(aes(fill=region, alpha = region), shape = 21, 
                size = 24, colour = "black", stroke = 0.8) +
  geom_tiplab(size=10, offset = 0.0000015)+
 
  theme_bw() +
  scale_fill_manual(values = c("grey80","royalblue4")) +
  scale_alpha_manual(values = c("Kenya-2025" = 1, "Global"=0.3))+
  
  xlab("Genetic distance") +
  ylab("Number of sequences") +
  labs(fill="Location") +
#  ylim(0,30) +
  theme(
    legend.position = "none", 
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.title = element_text(face = "bold", size = 40),
    legend.text = element_text(size = 28),
    axis.title.x = element_text(colour = "black", size = 40, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 36),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()) +
 # scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)), alpha="none") +
  geom_text2(aes(subset = !isTip, label = label),
             size = 9, color = "black", hjust = 1.5, vjust = -0.8) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.50)))



mpox_subset_plot

#ggsave("results/mpox_subset_plot.pdf", height = 45, width = 30)


################################################################################
#subset for 2024 samples 

mpox_subset.2024 <- tree_subset(mpox.ml.tree, "EPI_ISL_19302262/Kenya/2024-07-25",
                                levels_back=1)

mpox_subset_plot.2024 <- ggtree(mpox_subset.2024) %<+% metadata_final +
  theme_tree2() +
  geom_tippoint(aes(fill=region, alpha = region), shape = 21, 
                size = 24, colour = "black", stroke = 0.8) +
  geom_tiplab(size=10, offset = 0.0000015)+
  
  theme_bw() +
  scale_fill_manual(values = c("grey80","firebrick4")) +
  scale_alpha_manual(values = c("Kenya-2025" = 1, "Global"=0.3))+
  
  xlab("Genetic distance") +
  ylab("Number of sequences") +
  labs(fill="Location") +
  #  ylim(0,30) +
  theme(
    legend.position = "none", 
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.title = element_text(face = "bold", size = 40),
    legend.text = element_text(size = 28),
    axis.title.x = element_text(colour = "black", size = 40, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 36),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()) +
  # scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)), alpha="none") +
  geom_text2(aes(subset = !isTip, label = label),
             size = 9, color = "black", hjust = 1.5, vjust = -0.8) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.50)))



mpox_subset_plot.2024
#ggsave("results/mpox_subset_plot.2024.pdf", height = 45, width = 30)

##############################################################################
#  Final plot 
library(cowplot)

fig2 <- plot_grid(S1, mpox_plot, 
          mpox_subset_plot, 
          mpox_subset_plot.2024, scale = 0.95, 
          labels = c("A", "B", "C", "D"), 
          label_size = 40)


# fig2 <- ggdraw() +
#   draw_plot(S1, x=0, y=0.5, width = 0.48, height=0.45) +
#   draw_plot(mpox_plot, x = 0, y = 0, width = 0.48, height = 0.45) +
#   draw_plot(mpox_subset_plot, x = 0.52, y = 0.5, width = 0.48, height = 0.45) +
#   draw_plot(mpox_subset_plot.2024, x=0.52, y=0.0, width = 0.48, height = 0.45) +
#   draw_plot_label(label = c("A", "B", "C"), 
#                   x = c(0, 0.45, 0.45), 
#                   y = c(1, 0.95, 0.5), 
#                   size = 32, 
#                   fontface = "bold")


fig2
#ggsave("results/fig2.pdf", height = 38, width = 48, dpi = 300)
#gggsave("results/fig2.tiff", height = 38, width = 48, dpi = 300)



#################################################################################

