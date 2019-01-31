# Analyses for Figure 1
# January 2019
# Isla Myers-Smith

packrat::init()

# Packages ----
library(dplyr)
library(raster)
library(maps)
library(mapdata)
library(mapproj)
library(maptools)
library(rgdal)
library(spatial)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)

rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

# Load data ----
modisNDVIslope <- raster("data/MODIS_NDVI_slopes_2000-2015.tif")
gimmsNDVIslope <- raster("data/GIMMS_NDVI_slopes_2000-2015.tif")

# Reproject data
crslaea <- "+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

vegmap <- readOGR("shapefiles/cp_biozone_la_shp", "cavm_all polygon")
vegmap <- spTransform(vegmap, CRS(crslaea)) 
vegmap <- vegmap[which(vegmap$ZONE < 1),]

watermask <- raster("shapefiles/ocean_raster.tif")
water <- projectRaster(watermask, crs = crslaea)
water <- crop(water, extent(-3951586, 3945778, -3951599, 3950277))

arctic_circle <- readOGR("shapefiles/ArcticCircle", "ArcticCircle")
arctic_circle <- spTransform(arctic_circle, CRS(crslaea))

arctic_AE <- readOGR("shapefiles/ArcticAE", "ArcticAE_Combined_WGS")
arctic_AE <- spTransform(arctic_AE, CRS(crslaea))

tundra <- readOGR("shapefiles/Tundra65-78", "Tundra65-78")
tundra <- spTransform(tundra, CRS(crslaea))

EurAs <- readOGR("shapefiles/Eurasia180", "Eurasia180")
EurAs <- spTransform(EurAs, CRS(crslaea))

NorAm <- readOGR("shapefiles/NorthAmerica180", "NorthAmerica180")
NorAm <- spTransform(NorAm, CRS(crslaea))

# MODIS 2000 - 2015
modisNDVIslope <- projectRaster(modisNDVIslope, crs = crslaea)
modisNDVIslope <- (modisNDVIslope/10000)
modisNDVIslope.r <- modisNDVIslope
modisNDVIslope.r <- scale(modisNDVIslope.r)
r.min.modis <- cellStats(modisNDVIslope.r, "min")
r.max.modis <- cellStats(modisNDVIslope.r, "max")
modisNDVIslope.rs <- rescale(modisNDVIslope.r, x.min = r.min.modis, x.max = r.max.modis, new.min = -1, new.max = 1)

# GIMMS 2000 - 2015
gimmsNDVIslope <- projectRaster(gimmsNDVIslope, crs = crslaea)
gimmsNDVIslope.r <- gimmsNDVIslope
gimmsNDVIslope.r <- scale(gimmsNDVIslope.r)
r.min.gimms <- cellStats(gimmsNDVIslope.r, "min")
r.max.gimms <- cellStats(gimmsNDVIslope.r, "max")
gimmsNDVIslope.rs <- rescale(gimmsNDVIslope.r, x.min = r.min.gimms, x.max = r.max.gimms, new.min = -1, new.max = 1)

# Figures ----

# MODIS 2000 - 2015

tiff(file="plots/ArcticMODISNDVImap_2000-2015.tif", width=1000, height=1000)
par(mar=c(0,0,0,0), oma=c(6,6,6,6))

col5 <- colorRampPalette(c('#904E06', '#e2e2e2', '#2D821B'))
color_levels = 11
max_value = round(cellStats(modisNDVIslope, stat='max'),2)
color_sequence = round(seq(-1.5,1.5,length.out=color_levels+1),1)
col5(n=color_levels)

# Plot map
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(55, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
plot(laea_wrld_sp, col = 'black', lwd = 0.5)
plot(modisNDVIslope, add = T, col = col5(n=color_levels), breaks = c(-1.5, -0.05, -0.01, -0.005, -0.0025, -0.001, 0.001, 0.0025, 0.005, 0.01, 0.05, 1.5), axes=FALSE, legend = F)
plot(vegmap, add = T, col = 'white', border = FALSE)
plot(water, add = T, col= c("transparent", "white"), breaks = seq(0, 400, by = 200), legend = FALSE)
plot(arctic_circle, add = T, col= "black", lwd = 3, lty = "dashed")
plot(laea_wrld_sp, add = T, col = 'black', fill = FALSE, lwd = 0.5)
plot(arctic_AE, add = T, border='black', col = 'transparent', lwd = 3)
points(c(0,90), cex = 152, col = 'white', lwd = 160)
points(c(0,90), cex = 130, col = 'black', lwd = 1.5)

dev.off()

# GIMMS 2000 - 2015

tiff(file="plots/ArcticGIMMSNDVImap_2000-2015.tif", width=1000, height=1000)
par(mar=c(0,0,0,0), oma=c(6,6,6,6))

col5 <- colorRampPalette(c('#904E06', '#e2e2e2', '#2D821B'))
color_levels = 11
max_value = round(cellStats(gimmsNDVIslope, stat='max'),2)
color_sequence = round(seq(-1.5,1.5,length.out=color_levels+1),1)
col5(n=color_levels)

# Plot map
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(55, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
plot(laea_wrld_sp, col = 'black', lwd = 0.5)
plot(gimmsNDVIslope, add = T, col = col5(n=color_levels), breaks = c(-1.5, -0.05, -0.01, -0.005, -0.0025, -0.001, 0.001, 0.0025, 0.005, 0.01, 0.05, 1.5), axes=FALSE, legend = F)
plot(vegmap, add = T, col = 'white', border = FALSE)
plot(water, add = T, col= c("transparent", "white"), breaks = seq(0, 400, by = 200), legend = FALSE)
plot(arctic_circle, add = T, col= "black", lwd = 3, lty = "dashed")
plot(laea_wrld_sp, add = T, col = 'black', fill = FALSE, lwd = 0.5)
plot(arctic_AE, add = T, border='black', col = 'transparent', lwd = 3)
points(c(0,90), cex = 152, col = 'white', lwd = 160)
points(c(0,90), cex = 130, col = 'black', lwd = 1.5)

dev.off()

tiff(file="plots/world_map.tif", width=1000, height=1000)
par(mar=c(0,0,0,0), oma=c(6,6,6,6))

# Plot map
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(55, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
plot(laea_wrld_sp, col = 'black', lwd = 0.5)
plot(vegmap, add = T, col = 'white', border = FALSE)
plot(water, add = T, col= c("transparent", "white"), breaks = seq(0, 400, by = 200), legend = FALSE)
plot(arctic_circle, add = T, col= "black", lwd = 3, lty = "dashed")
plot(laea_wrld_sp, add = T, col = 'black', fill = FALSE, lwd = 0.5)
plot(tundra, add = T, border='black', col = '#00705280', lwd = 3)
plot(EurAs, add = T, border='black', col = '#e8a60080', lwd = 3)
plot(NorAm, add = T, border='black', col = '#4286f480', lwd = 3)
points(c(0,90), cex = 152, col = 'white', lwd = 160)
points(c(0,90), cex = 130, col = 'black', lwd = 1.5)

dev.off()
