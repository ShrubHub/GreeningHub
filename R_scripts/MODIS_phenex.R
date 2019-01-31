# Greenup function for Greening Note
# December 2018
# Isla Myers-Smith and Jeff Kerby

packrat::init()

# Libraries ----
library(phenex)
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)

# DOY Function ----

DOY.data <- function(year, DOY, NDVI) {
  y <- cbind.data.frame(DOY, NDVI)
  names(y) <- c("DOY", "NDVI")
  if(leapYears(year[1]) == TRUE){
    x <- cbind.data.frame(seq(1, 366))
    names(x) <- c("DOY")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  } else {
    x <- cbind.data.frame(seq(1, 365))
    names(x) <- c("DOY")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  }
  return(z)
}

# Read data ----
MODIS <- read.csv("data/modisData_QHI_Kanger.csv")
MODIS <- MODIS %>% mutate(year = gsub("-.*$", "", as.character(MODIS$date)), DOY = as.integer(strftime(date, format = "%j")))
MODIS$year <- as.numeric(as.character(MODIS$year))
MODIS$SummaryQA <- as.numeric(as.character(MODIS$SummaryQA))
MODIS$NDVI <- as.numeric(as.character(MODIS$NDVI))
MODIS$EVI <- as.numeric(as.character(MODIS$EVI))
MODIS$NDVI <- as.numeric(as.character(MODIS$NDVI))

MODIS <- MODIS %>% filter(!is.na(NDVI) | SummaryQA <= 3) %>% mutate(NDVI=NDVI/10000, EVI=EVI/10000)

# Reformat MODIS data
MODIS.data <- cbind.data.frame(site = as.character(MODIS$site), NDVI = as.numeric(MODIS$NDVI), DOY = MODIS$DOY, year = MODIS$year) 

# Remove negative values
# MODIS.data <- MODIS.data %>% filter(MODIS.data$NDVI > 0)

# Add zeros
num.years <- max(MODIS.data$year) - min(MODIS.data$year) + 1
num.sites <- length(unique(MODIS.data$site))

# Phenex calculations ----

DOY <- cbind.data.frame(DOY = rep(c(rep(1, num.years), rep(32, num.years), rep(60, num.years), rep(305, num.years), rep(335, num.years)), num.sites))
NDVI <- cbind.data.frame(NDVI = rep(0, length(DOY$DOY)))
year <- cbind.data.frame(year = rep(rep(seq(min(MODIS.data$year),max(MODIS.data$year)), 5), num.sites))
site <- cbind.data.frame(site = rep(rep(as.character(unique(MODIS.data$site)), 5), num.years))
site <- arrange(site, site)
zero.data <- cbind.data.frame(site, year, DOY, NDVI)

MODIS.data <- rbind(MODIS.data, zero.data)

# run DOY function to create vectors for modelNDVI
greenup <- MODIS.data %>% group_by(site, year) %>% do(., DOY.data(.$year, .$DOY, .$NDVI))

# run modelNDVI function - takes a while
greenup.all <- greenup %>% 
  # group by site and year
  group_by(site, year) %>% 
  # apply modelNDVI function
  do(., ndvi.values = unlist(modelNDVI(ndvi.values=.$NDVI, year.int=.$year, correction='none', method="DLogistic", MARGIN=2, multipleSeasons=FALSE, doParallel=FALSE, silent=TRUE))) %>%
  # Save just the DOY and modelledValues in individual columns of the tibble for plotting
  mutate(modVal = ifelse(leapYears(year[1])==TRUE, 
                         list(cbind(DOY = seq(1:366), modVal = modelledValues(ndvi.values[[1]]))), 
                         list(cbind(DOY = seq(1:365), modVal = modelledValues(ndvi.values[[1]])))
  )) %>% 
  # create phenoPhase dates and max NDVI
  mutate(
    greenup.date.05 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.05)$mean, 
    greenup.date.50 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.50)$mean, 
    greenup.date.95 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.95)$mean, 
    ndvi.date.max = phenoPhase(ndvi.values[[1]], phase="maxval", method="local")$mean, 
    senescence.date.05 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.05)$mean, 
    senescence.date.50 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.50)$mean, 
    senescence.date.95 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.95)$mean,
    gs.length.05 = senescence.date.05 - greenup.date.05,
    gs.length.50 = senescence.date.50 - greenup.date.50,
    gs.length.95 = senescence.date.95 - greenup.date.95
    #, integrateTimeserie = integrateTimeserie(ndvi.values[[1]][[1]], start=greenup.date.50, end=senescence.date.50, n=1000)
  )

# Save the greenup df for faster loading
save(greenup.all, file="data/greenup.qhi.kanger.RData")
