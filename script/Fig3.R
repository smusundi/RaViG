
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


loadfonts(device = "pdf")

# colors assigned 
cols <- c("gold", "firebrick3", "orange",  "aquamarine4",  "blueviolet",
          "salmon", "olivedrab3","royalblue4",  "lightpink4", "black",
          "powderblue",  
          "cadetblue1", "darkolivegreen1", "deeppink1", "gray")



utange_chickenpox <-  read.tree("data/output.chickenpox.treefile")


utange_chickenpox$tip.label

# GenBank metadata  for VZVZ 
metadata_chickenpox <- read_tsv("data/chickenpox_metadata.tsv") %>%
  filter(Length>100000) %>%
  select(1,6, `Isolate Collection date`) %>%
  na.omit()

colnames(metadata_chickenpox)[1:3] <- c("accession", "location", "date")

str(metadata_chickenpox)




kenya_vzv <- data.frame(accession=c("RVG-005", "RVG-007", "RVG-020", "RVG-032"), 
                        location= "Kenya", 
                        date=c("2025-11-19","2025-11-19", "2025-11-20", 
                               "2025-11-25"))

genebank.vzv.metadata <- rbind(metadata_chickenpox, kenya_vzv)


#################################################################################
# Nextclade metadata 
nextclade_VZV.metadata <- read_tsv("data/VZV_nextclade.tsv") %>%
  select(seqName, clade)


colnames(nextclade_VZV.metadata) [1:2] <- c("accession", "clade")


kenya_clades <- data.frame(accession=c("RVG-020", "RVG-032", "RVG-005", "RVG-007"), 
                           clade= "clade 5")

nextclade.vzv <-  rbind(nextclade_VZV.metadata, kenya_clades)


combined.vzv.metadata <- genebank.vzv.metadata %>%
  left_join(nextclade.vzv, by="accession")



# drop tips lacking metadata 
dropped_tips <- setdiff(utange_chickenpox$tip.label, combined.vzv.metadata$accession)



# state outliers / remove tips with very long branches 

long_outliers <- c("PP378487.1", "PP378488.1", "PP378489.1", "MH379685.1", 
                   "KU926318.1", "KU926317.1", "MT370830.1", "KU926320.1", 
                   "PV061461.1", "MH499468.1", "OQ835722.1")


tips.dropped <- c(dropped_tips, long_outliers)



utange_chickenpox <-  ape::drop.tip(utange_chickenpox, tips.dropped)



length(utange_chickenpox$tip.label)
############################################################################
# modify names of the tree files to include country and dates for sample collection

tree_ids <- combined.vzv.metadata %>%
  filter(accession %in% utange_chickenpox$tip.label) %>%
  mutate(sampleid=paste(accession, location, date, sep = "/")) 


new_labels <- tree_ids$sampleid[match(utange_chickenpox$tip.label, tree_ids$accession)]

utange_chickenpox$tip.label <-  ifelse(is.na(new_labels), utange_chickenpox$tip.label, 
                                       new_labels)


utange_chickenpox$tip.label

###################################################################

utange_chickenpox <- midpoint(utange_chickenpox)



# filter metadata 
combined.vzv.metadata <- tree_ids %>%
  filter(sampleid %in% utange_chickenpox$tip.label) %>%
  mutate(Location=case_when(str_detect(accession, "RVG") ~ "Kenya-2025",
                            str_detect(accession, "PQ505474.1") ~ "Kenya-2024", 
                            location %in% c("Germany", "United Kingdom", "Italy", "Finland", 
                                            "Spain", "Sweden", "Russia", "Portugal", "France", 
                                            "Greece") | str_detect(location, "Germany") ~ "Europe", 
                            location %in% c("India", "South Korea", "China", "Pakistan", "Singapore") | str_detect(location, "India") ~ "Asia", 
                            location %in% c("USA", "Mexico") | str_detect(location, "USA")~ "Americas", 
                            .default = "Africa"), 
         Location=as.factor(Location), 
         region=case_when(Location == "Kenya-2025" ~ "Kenya-2025", 
                          Location == "Kenya-2024" ~ "Kenya-2024", 
                          .default = "Global"), 
         clade=as.factor(clade), location=as.factor(location), 
         region=as.factor(region)) %>%
  select(sampleid, accession, location, clade, date, Location, region)


#write.tree(file = "data/tree.nwk", utange_chickenpox)
# Check number of sequences deposited per country across time

S2 <- combined.vzv.metadata %>%
  mutate(date = case_when( nchar(date) == 4 ~ paste0(date, "-01-01"),
                           nchar(date) == 7 ~ paste0(date, "-01"),
                           TRUE ~ date),date = as.Date(date),
         month_date=floor_date(date, unit = "month"),
         location=str_remove_all(location, ":.*$"),
         location=as.factor(location)) %>%
  select(sampleid, location, month_date) %>%
  group_by(location, month_date) %>%
  summarise(count=n(), .groups = "drop") %>%
  mutate(country = fct_reorder(location, month_date, .fun = min)) %>%
  ggplot(aes(month_date, location, size = count, color = count)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "Sequences", range = c(2, 30)) +
  scale_color_viridis_c(name = "Sequences") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "5 years") +
  labs(x = "Time", y = "Country") +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = "black", size = 24),
    axis.text.x = element_text( angle=90, hjust = 1, colour = "black", size = 24),
    axis.title = element_text(colour = "black", face = "bold", size = 30), 
    legend.title = element_text(size = 30, face = "bold", vjust = 1.5), 
    legend.text = element_text(size = 24), 
    legend.spacing.y = unit(5, "cm")) +
  guides(colour = guide_colourbar(barheight = unit(10, "cm"),
                                  barwidth  = unit(1, "cm")))




# plot divergence tree 
Fig3 <-  ggtree(utange_chickenpox)  %<+%  combined.vzv.metadata +
  theme_tree() +
  geom_tippoint(aes(fill = region, alpha = region), size=16.0,  colour = "black", shape = 21,
                stroke = 0.8) +
  #scale_shape_manual(values = 15:24)+
  #geom_text(aes(label = node), size=2) +
  # geom_tiplab(size=2)+
  #geom_treescale(y = -5, offset = -1, fontsize = 3, linesize = 0.4) +
  #geom_tiplab(aes(node = node, 
  #                color = ifelse(grepl("RVG", label), "Kenyan", "Global")),
  #            geom = "text", fontface = 2.5, size = 2.5, show.legend = FALSE) +
  #scale_fill_manual(values = cols) +
  scale_fill_manual(values = c("grey80", "firebrick", "royalblue4"))+
  scale_alpha_manual(values = c("Kenya-2025" = 1, "Kenya-2024" = 1, "Global" = 0.6))+
  # labs(color = "Location") +
  guides(fill = guide_legend(override.aes = list(alpha = 1)), alpha="none") +
  xlab("Genetic distance") +
  ylab("Number of sequences") +
  labs(fill="Region") +
  theme_bw() +
  ylim(0,265) +
  theme( legend.position = c(0.08, 0.92),
         legend.justification = c(0, 1),
         legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
         legend.title = element_text(face = "bold", size = 24),
         legend.text = element_text(size = 20),
         axis.title.x = element_text(colour = "black", size = 24, face = "bold"),
         axis.text.x = element_text(colour = "black", size = 20),
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         panel.grid = element_blank()) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  geom_cladelab(node = 306, label = "Clade 5", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.2, align = TRUE) +
  geom_cladelab(node = 451, label = "Clade 3", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 312, label = "Clade 1", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 441, label = "Clade 9", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 476, label = "Clade 6", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 347, label = "Clade 2", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 376, label = "Clade 2 vaccine", color = "black", offset = .00001,
                fontsize = 6.5, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = c(415, 349), label = "Clade 2 wild type", color = "black", offset = .00001,
                fontsize = 6.5, barsize = 1.5, align = TRUE) +
  geom_cladelab(node = 427, label = "Clade 4", color = "black", offset = .0005,
                fontsize = 8, barsize = 1.5, align = TRUE)


Fig3
###############################################################################
# Subplot for VZV 
vzv.subset.005 <-  tree_subset(utange_chickenpox,  "PQ505474.1/Kenya/2024-08-23",
                               levels_back =20)


rvg005.plot <- ggtree(vzv.subset.005)  %<+%  combined.vzv.metadata +
  theme_tree() +
  geom_tippoint(aes(fill = region, alpha = region), size=12.0, shape = 21, colour = "black", 
                stroke = 0.8)+
  geom_tiplab(size=6, offset = 0.00001) +
  scale_fill_manual(values = c("grey80", "firebrick", "royalblue4"))+
  scale_alpha_manual(values = c("Kenya-2025" = 1, "Kenya-2024" = 1, "Global" = 0.6))+
  labs(color = "Location") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 12))) +
  xlab("Genetic distance") +
  ylab("Number of sequences") +
  labs(fill="Region") +
  theme_bw() +
  theme( axis.title.x = element_text(colour = "black", size = 24, face = "bold"),
         axis.text.x = element_text(colour = "black", size = 20),
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         panel.grid = element_blank(), 
         legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.5))) +
  # scale_y_continuous(expand = expansion(mult = c(0.01, 0.009))) +
  geom_text2(aes(subset = !isTip, label = label),
             size = 4.5, color = "black", hjust = 1.5, vjust = -0.8)

rvg005.plot

Fig3A <- plot_grid(S2, Fig3, 
          nrow = 2, labels = c("A", "B"), label_size = 40, scale = 0.95)

Fig3B <-  plot_grid(Fig3A, rvg005.plot, labels = c("", "C"), 
                    label_size = 40, scale = 0.95)

# 
# 
# fig3<- ggdraw() +
#   draw_plot(Fig3, x = 0, y = 0, width = 0.45, height = 1) +
#   draw_plot(rvg005.plot, x = 0.5, y = 0, width = 0.5, height = 1) +
#   draw_plot_label(label = c("A", "B"), 
#                   x = c(0, 0.5), 
#                   y = c(1, 1), 
#                   size = 32, 
#                   fontface = "bold")




#ggsave("results/Fig3B.tiff", height = 36, width = 32, dpi = 400)

##############################################################################

  
