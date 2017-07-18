require(ggplot2)
require(ggmap)
require(maps)
require(mapproj)
require(mapdata)
require(rgeos)
require(maptools)
require(sp)
require(raster)
require(rgdal)
require(dismo)

source("analysis/InitialSetup.R")

range(env$latitude)
range(env$longitude)

latrange <- c()

basemap <- get_map(location = c(-122.3, 44.2, -122.1, 44.3), 
                   zoom = 12, maptype = "terrain")
map1 <- ggmap(basemap)
map1 + geom_point(data = env, 
                  aes(x = longitude, y = latitude, cex = 2), 
                  color = "red", show.legend = F) +
  geom_polygon()
  scale_x_continuous(limits = c(-122.27, -122.10)) +
  scale_y_continuous(limits = c(44.18, 44.30))

streamlines <- readShapeSpatial("imagery/lidar_stream/lidar_stream.shp", 
                 proj4string = CRS("+proj=longlat +datum=WGS84"))
streamlines <- project(streamlines, "+proj=longlat")
data <- fortify(streamlines)
ggmap(map1)
