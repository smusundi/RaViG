
# clear working environment
rm(list=ls())

# load packages 
# If the following packages are not installed please use the install.library("package name")
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
cols <- c("gold4", "firebrick3", "orange4",  "aquamarine3",  "blueviolet",
          "salmon3", "olivedrab3","royalblue", "wheat2", 
          "powderblue", "lightpink4", "coral", 
          "cadetblue1", "darkolivegreen", "deeppink1", "gray")

###############################################################################
# Generate maps -
################################################################################

# ── 1. LOAD & PREP SUBCOUNTY DATA ─────────────────────────────────────────────

mbaklf <- st_read("data/geoBoundaries-KEN-ADM2-all/geoBoundaries-KEN-ADM2_simplified.shp")

coastal_subcounties <- c("Kisauni", "Likoni", "Jomvu", "Changamwe",
                         "Mvita", "Nyali", "Rabai", "Kaloleni")

mbaklf <- mbaklf %>%
  filter(shapeName %in% coastal_subcounties)

study_counts <- data.frame(
  shapeName    = c("Kisauni", "Likoni", "Changamwe", "Kaloleni", "Jomvu", "Mvita"),
  participants = c(23, 3, 2, 2, 1, 1))

mbaklf_map <- mbaklf %>%
  left_join(study_counts, by = "shapeName")


# ── 2. LOAD KENYA COUNTIES ────────────────────────────────────────────────────

# Use ADM1 from geoBoundaries for consistency with your data source
kenya_adm1 <- st_read("data/geoBoundaries-KEN-ADM1-all/geoBoundaries-KEN-ADM1_simplified.shp")

# check names
unique(kenya_adm1$shapeName)

# tag coastal county
kenya_adm1 <- kenya_adm1 %>%
  mutate(fill_group = case_when(
    shapeName == "Mombasa" ~ "Mombasa",
    shapeName %in% c("Kilifi", "Kwale", "Taita-Taveta",
                     "Lamu", "Tana River") ~ "Coastal",
    TRUE ~ "Other"))


# ── 3. KENYA INSET ────────────────────────────────────────────────────────────

# bounding box of your subcounty map (to draw zoom box on inset)
bbox <- st_bbox(mbaklf_map)

kenya_inset <- ggplot() +
  
  geom_sf(data = kenya_adm1,
          aes(fill = fill_group),
          color = "grey50", linewidth = 0.2) +
  
  scale_fill_manual(
    values = c(
      "Mombasa" = "grey70",     # dark purple — matches viridis
      "Coastal" = "grey70",     # light green — coastal context
      "Other"   = "white"),
    guide = "none") +
  
  # zoom box showing extent of main map
  annotate("rect",
           xmin = bbox["xmin"] - 0.02,
           xmax = bbox["xmax"] + 0.02,
           ymin = bbox["ymin"] - 0.02,
           ymax = bbox["ymax"] + 0.02,
           fill = NA, color = "black", linewidth = 0.8) +
  
  coord_sf() +
  theme_void() +
  theme( panel.background = element_rect(fill = "grey97", color = "black",
                                         linewidth = 0.7),
         plot.background  = element_rect(fill = "grey97", color = "black",
                                         linewidth = 0.7))


# ── 4. MAIN SUBCOUNTY MAP ─────────────────────────────────────────────────────

mbaklf_plot <- ggplot(mbaklf_map) +
  
  geom_sf(aes(fill = participants),
          color = "grey25", linewidth = 0.4) +
  
  scale_fill_viridis_c(
    option = "viridis", direction = -1,
    na.value = "grey95", name = "Participants",
    breaks = c(5, 10, 15, 20),
    guide = guide_colorbar(
      barwidth = 1, barheight = 6,
      title.position = "top", title.hjust = 0.5)) +
  
  geom_text_repel(
    data = mbaklf_map,
    aes(label = shapeName, geometry = geometry),
    stat = "sf_coordinates",
    size = 3.2, fontface = "bold", color = "black",
    box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey40", segment.size = 0.3,
    min.segment.length = 0,
    bg.color = "white", bg.r = 0.15) +
  
  annotation_scale(
    location = "bl", width_hint = 0.25,
    text_cex = 0.7, line_width = 0.5) +
  
  annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(1, "cm"), width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(text_size = 8)) +
  
  coord_sf() +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.position = "right",   # bottom left — away from inset
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8))

mbaklf_plot
# ── 5. COMBINE ────────────────────────────────────────────────────────────────

final_map <- ggdraw(mbaklf_plot) +
  
  draw_plot(
    kenya_inset,
    x      = 0.01,   # top right corner
    y      = 0.80,
    width  = 0.10,
    height = 0.10) 

final_map



############################################################################
# Ct values and number of reads 
options(scipen = 1000)

ravig_data_reads <-  read_excel("data/ravig_data_reads.xlsx", sheet = 2)

colnames(ravig_data_reads)

VZV_ct <- ravig_data_reads %>%
  mutate(MPXV_Ct_values=as.numeric(MPXV_Ct_values), 
         average=as.numeric(average), 
         vzv_infection=case_when(!is.na(MPXV_Ct_values) & average>=1 ~ "VZV + MPXV", 
                                 is.na(MPXV_Ct_values) & average>=1 ~ "VZV", 
                                 MPXV_Ct_values>1 & is.na(average) ~ "MPXV", 
                                 is.na(MPXV_Ct_values) & is.na(average) ~ "Neg")) %>%
  select(vzv_infection,average) %>%
  na.omit() %>%
  suppressWarnings()


ct_plot <- ggplot(VZV_ct , aes(x = vzv_infection, y = average)) +
  geom_boxplot(size = 0.2, width=0.5) +
  geom_jitter(shape=21, position=position_jitter(0.08), size=3, fill="gold3") +
  ylab("Ct values") +
  xlab("VZV Infection") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"), 
        axis.title = element_text(colour = "black", face = "bold"))

ct_plot

ct_vs_reads <- ggplot(ravig_data_reads, aes(average, VZV_reads)) +
  geom_point(size=3, fill="gold3", shape=21) +
  xlab("Ct values") +
  ylab("VZV reads (log10 scale)") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"), 
        axis.title = element_text(colour = "black", face = "bold")) +
  scale_y_log10()

ct_vs_reads

ct_values_coverage <- plot_grid(ct_plot, ct_vs_reads, 
                                labels = c("B", "C"), scale = 0.75, nrow = 2)

Fig1<- plot_grid(final_map,ct_values_coverage, 
                 labels = c("A", ""))

Fig1



#ggsave("results/Fig1.tiff", height = 8, width = 10, dpi = 300)

