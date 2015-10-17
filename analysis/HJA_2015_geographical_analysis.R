# Load packages
require("vegan")
require("sp")
require("geoR")
require("rgdal")
require("raster")
require("RgoogleMaps")
require("maptools")

setwd("~/GitHub/HJA-streams/")

# Load the environmental data
hja <- read.table(file="./data/water_chemistry//2015-07-27_HJA-env.csv", head=TRUE, sep=",")

# Parse out latitudes and longitudes
lat <- as.numeric(hja[, 4])
long <- as.numeric(hja[, 5])

# Map samples and data for LC
HJAmap <- GetMap(center = c(lat=44.24,lon=-122.17), zoom = 12,
                 destfile = "HJAmap.png", maptype = "terrain")

PlotOnStaticMap(HJAmap, zoom = 12, cex=2, col='blue')
PlotOnStaticMap(HJAmap, lat, long, cex=1, pch=20, col='red', add=TRUE)

WS01map <- GetMap(center = c(lat=44.204, lon=-122.25), zoom = 15,
                  destfile = "WS01map.png", maptype = "terrain")
PlotOnStaticMap(WS01map, zoom = 12, cex = 2, col = 'blue')
PlotOnStaticMap(WS01map, lat, long, cex=1, pch=20, col='red', add=TRUE)

# Plot Variogram
long_east = 360 + hja$Longitude

xy <- data.frame(name = hja$Name, lat = hja$Latitude, long = long_east)
coordinates(xy) <- c("long", "lat")
proj4string(xy) <- CRS("+proj=longlat +datum=NAD83")
UTM <- spTransform(xy, CRS("+proj=utm +zone=10 +ellps=WGS84"))

UTM <- as.data.frame(UTM)
coordsUTM <- UTM[,2:3]

dev.off()
var.trend.geoR <- variog(coords=coordsUTM, data=lc$Elevation, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for Elevation (m)",
     xlab = "Distance (m)", col="blue")

var.trend.geoR <- variog(coords=coordsUTM, data=lc$SpCond, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for SpC",
     xlab = "Distance (m)", col="blue")

var.trend.geoR <- variog(coords=coordsUTM, data=lc$pH, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for pH",
     xlab = "Distance (m)", col="blue")

var.trend.geoR <- variog(coords=coordsUTM, data=lc$ORP, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for ORP",
     xlab = "Distance (m)", col="blue")

var.trend.geoR <- variog(coords=coordsUTM, data=lc$TDS, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for TDS",
     xlab = "Distance (m)", col="blue")

var.trend.geoR <- variog(coords=coordsUTM, data=lc$Temperature, option="bin")
plot(var.trend.geoR, type = "b", main = "Variogram for Temp",
     xlab = "Distance (m)", col="blue")

# Layers
LandCover <- raster("~/Box Sync/2015_Andrews/Imagery/LandCover_2011/NLCD_OR_20111.ovr")
plot(LandCover, xlab="Longitude", ylab="Latitude",
     main="Map of landcover")

Geology <- readShapeSpatial("~/Box Sync/2015_Andrews/Imagery/geo_gis/DOGAMI//OGDCv5//shp//G_REF_MAP.shp")
plot(Geology)
Geology.unit <- readShapeSpatial("~/Box Sync/2015_Andrews/Imagery/geo_gis/DOGAMI//OGDCv5//shp//G_MAP_UNIT.shp")
plot(Geology.unit)

climate <- read.table("./data/hja_site_data/climate.dat", sep='\t',
                      na.strings="")
