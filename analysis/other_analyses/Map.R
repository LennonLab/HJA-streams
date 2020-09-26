source("analysis/InitialSetup.R")

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

range(env$latitude)
range(env$longitude)

latrange <- c()

basemap <- get_map(location = c(-122.3, 44.2, -122.1, 44.3), 
                   zoom = 12, maptype = "terrain")
map1 <- ggmap(basemap)
hja.map <- map1 + geom_point(data = env, 
                  aes(x = longitude, y = latitude, cex = 2), 
                  color = "black", shape = 21, fill = "red", show.legend = F) +
  scale_x_continuous(limits = c(-122.27, -122.10)) +
  scale_y_continuous(limits = c(44.18, 44.30))+
  labs(x="Longitude (dec. degrees)", y="Latitude (dec. degrees)")+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
hja.map
