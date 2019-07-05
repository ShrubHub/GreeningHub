# Analyses for Figure 1
# January 2019
# Isla Myers-Smith

#packrat::init()

# Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)

# Data ----

# Phenology
bands <- read.csv("data/NDVI_bands.csv", stringsAsFactors = F)

pallet <- viridis_pal()(20)
pallet2 <- c(pallet[20], pallet[18], pallet[17], pallet[16], pallet[15], pallet[10], pallet[6], pallet[1])
pallet3 <- c("#ca0000", "#faff4a", "#ffc900", "#ff9600", "#ff6000", "#009933", "#02bef1", "#4600aa")
bands$Satellite <- factor(bands$Satellite, levels = c("AVHRR", "Landsat1-5-MSS", "Landsat1-5-TM", "Landsat7-ETM+", "Landsat8-OLI-TIRS", "MODIS", "Sentinel-2A", "Sentinel-2B"))

# Figures ----

plot <- ggplot(bands) + 
  geom_rect(aes(xmin = min_band, xmax = max_band, ymin = min_year, ymax = max_year, colour = Satellite, fill = Satellite), alpha = 0.1, lwd=0.5) +
   ylab("Deployment (years)\n") +
   xlab("Wavelength (nm)") +
   scale_x_continuous(breaks = c(500, 600, 700, 800, 900, 1000, 1100, 1200), limits = c(500, 1200)) +
   scale_fill_manual(values = pallet3) +
   scale_colour_manual(values = pallet3) +
   theme_bw() +
   theme(axis.text = element_text(size = 14), 
         axis.title = element_text(size = 16),
         axis.text.x = element_text(angle = -45, hjust = -0.05),
         axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),                                          
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_blank(),  
         plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
         plot.title = element_text(size=16, vjust=1, hjust=0),
         legend.text = element_text(size=14),          
         legend.title = element_blank(),                              
         legend.key = element_blank(),
         legend.spacing.x = unit(0.5, 'cm'),
         legend.key.size = unit(1.5, 'lines'),
         legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))

ggsave("plots/band_widths.png", plot = plot, width = 10, height = 4, units = "in", dpi = 300)
