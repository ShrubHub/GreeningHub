# Analyses for Figure 3
# January 2019
# Isla Myers-Smith

packrat::init()

# Packages ----
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)

# Load data ----

sat_data <- read.csv("data/Qikiqtaruk_vegtype_satellite_data.csv")
sat_data <- dplyr::select(sat_data, -time, -id, -longitude, -latitude)
sat_data_mean <- sat_data %>% group_by(data, year, vegType, date) %>% summarise(NDVI = mean(NDVI))

drone_data <- read.csv("data/Qikiqtaruk_vegtype_drone_data.csv")
drone_data_mean <- drone_data %>% group_by(data, year, vegType, date) %>% summarise(NDVI = mean(NDVI))

data <- full_join(sat_data, drone_data_mean)
data <- distinct(data)
data$date <- data$date %>% str_replace("date*", "")
data$data <- factor(data$data, levels = c("Drone", "Sentinel", "Landsat", "MODIS"))
data <- data %>% group_by(data, year, date, vegType, NDVI) %>% mutate(datayear = paste0(data, " ", year))
data$datayear <- factor(data$datayear, levels = c("Drone 2017", "Sentinel 2016", "Sentinel 2017", "Landsat 2016", "Landsat 2017", "MODIS 2016", "MODIS 2017"))
data_mean <- data %>% group_by(data, year, vegType, datayear) %>% summarise(NDVI = mean(NDVI)) 

# Figures ----

all <- ggplot(data) +
  geom_point(aes(x = datayear, y = NDVI, colour = vegType), alpha = 0.8, size = 2, show.legend = F) +
  ylab("NDVI\n") +
  xlab("") +
  scale_colour_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  scale_fill_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  ylim(0.3,0.9) +
  annotate("text", x = 0.6, y = 0.9, label = "A. NDVI data across scales", size=5, hjust=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks.length=unit(.25, "cm"))

all <- all + geom_label_repel(
  aes(datayear, NDVI, label = date), size = 2,
  box.padding = 0.3, point.padding = 0.5, 
  segment.color = 'grey50')

drone <- ggplot(drone_data) +
  geom_violin(aes(x = vegType, y = NDVI, colour = vegType, fill = vegType), alpha = 0.8, size = 0.5) +
  ylab("NDVI\n") +
  xlab("") +
  scale_colour_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  scale_fill_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  ylim(0.3,0.9) +
  annotate("text", x = 0.6, y = 0.9, label = "B. Drone NDVI data", size=5, hjust=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks.length=unit(.25, "cm"))

plot <- grid.arrange(all, drone, ncol = 2, widths=c(8,4))

ggsave("plots/NDVI_GEE_vegType_detailed.png", plot = plot, width = 12, height = 6, units = "in", dpi = 300)

data_mean_2017 <- data_mean %>% filter(year == 2017) %>% mutate(vegType2 = vegType)
levels(data_mean_2017$vegType2) <- c(levels(data_mean_2017$vegType2), "Tussock Tundra", "Grass/Herb Tundra", "Shrub Tundra") 
data_mean_2017$vegType2[data_mean_2017$vegType2 == "Herschel"] <- "Tussock Tundra"
data_mean_2017$vegType2[data_mean_2017$vegType2 == "Komakuk"] <- "Grass/Herb Tundra"
data_mean_2017$vegType2[data_mean_2017$vegType2 == "Shrub"] <- "Shrub Tundra"

mean <- ggplot(data_mean_2017) +
  geom_point(aes(x = datayear, y = NDVI, colour = vegType2, fill = vegType2, size = data), alpha = 0.5, shape = 22) +
  ylab("NDVI\n") +
  xlab("") +
  scale_colour_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  scale_fill_manual(values = c("#4286f4", "#e8a600", "#007052")) +
  scale_size_manual(values = c(3, 6, 9, 14)) +
  ylim(0,1) +
  #annotate("text", x = 0.6, y = 1, label = "A. Peak growing season NDVI", size=5, hjust=0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks.length=unit(.25, "cm"))

ggsave("plots/NDVI_vegType.png", plot = mean, width = 5, height = 4, units = "in", dpi = 300)
