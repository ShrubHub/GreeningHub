# Analyses for Figure 1
# January 2019
# Isla Myers-Smith

packrat::init()

# Packages ----
library(tidyverse)
library(ggplot2)
library(gridExtra)

# Load data ----
modis <- read.csv("data/MODIS_NDVItimeseries_QAs_regions.csv")
modis$NDVI <- modis$NDVI/10000
modis <- na.omit(modis)

gimms <- read.csv("data/GIMMS_NDVItimeseries_regions.csv")

modis$region <- factor(modis$region, levels=c("Eurasia", "Tundra", "North America"))
gimms$region <- factor(gimms$region, levels=c("Eurasia", "Tundra", "North America"))

# Figures ----

modis.plot <- ggplot(modis) +
  geom_point(aes(x = year, y = NDVI, colour = region), alpha = 0.5, size = 2) +
  geom_line(aes(x = year, y = NDVI, colour = region), alpha = 0.5, size = 0.5) +
  geom_smooth(method=lm, aes(x = year, y = NDVI, colour = region, fill = region), alpha = 0.2, show.legend=F) + 
  ylab("Annual Max. NDVI\n") +
  xlab("") +
  #xlim(1982, 2016) +
  scale_x_continuous(limits = c(1982, 2018), breaks = c(1982, 1990, 2000, 2010, 2018)) +
  ylim(0.3, 0.8) +
  scale_colour_manual(values = c("#e8a600", "#007052", "#4286f4")) +
  scale_fill_manual(values = c("#e8a600", "#007052", "#4286f4")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))

gimms.plot <- ggplot(gimms) +
  geom_point(aes(x = year, y = NDVI, colour = region), alpha = 0.5, size = 2) +
  geom_line(aes(x = year, y = NDVI, colour = region), alpha = 0.5, size = 0.5) +
  geom_smooth(method=lm, aes(x = year, y = NDVI, colour = region, fill = region), alpha = 0.2, show.legend=F) + 
  ylab("Annual Max. NDVI\n") +
  xlab("") +
  #xlim(1982, 2016) +
  scale_x_continuous(limits = c(1982, 2018), breaks = c(1982, 1990, 2000, 2010, 2018)) +
  ylim(0.3, 0.8) +
  scale_colour_manual(values = c("#e8a600", "#007052", "#4286f4")) +
  scale_fill_manual(values = c("#e8a600", "#007052", "#4286f4")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))

plot <- grid.arrange(gimms.plot, modis.plot, nrow = 2)

ggsave("plots/NDVI_trends.png", plot = plot, width = 7, height = 7, units = "in", dpi = 300)
