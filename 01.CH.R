options(java.parameters = "-Xmx6g")

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(sf)
library(maptools)
library(dplyr)
library(fasterize)

wdpa_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"


# Reading PA dataset
pa <- raster("PA_1km_up.tif") 
crs(pa) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"


# Land Map
shp <- raster("LandMap_resampled_up.tif")
crs(shp) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

# Reading species dataset
dc_cl <- read.csv("Data_compiled.csv", header = T)

dc_cl <- dc_cl %>%
  dplyr::select("species", "decimalLon", "decimalLat")


# Remove records without coordinates
dc_cl <- dc_cl%>%
  filter(!is.na(decimalLon))%>%
  filter(!is.na(decimalLat))%>%
  filter(!is.na(species))


# Remove duplicated records
dc_cl <- dc_cl[!duplicated(dc_cl),]

# Remove species with low occurrence records
dc_cl <- dc_cl %>%
  group_by(species) %>%
  filter(n() > 2) %>%
  ungroup()

species <- unique(dc_cl$species)


for (i in species) {
  print(i)
  speciesname <- gsub(" ", "_", i)
  dc <- dc_cl %>% filter(species==i)
  dc <- dc %>%
    dplyr::select("decimalLon", "decimalLat")
  
  # Convex Hull
  coords <- na.omit(dc[,c('decimalLon','decimalLat')])
  dc <- dc[!duplicated(paste(dc$decimalLat, dc$decimalLat)), ]
  coordinates(coords) <- ~decimalLon + decimalLat
  hull <- gConvexHull(coords)
  hull <- as(hull, "sf")
  hull <- st_make_valid(hull)
  hull <- st_set_crs(hull, 4326)
  hull <- st_transform(hull, crs = wdpa_crs)
  
  
  # Rasterize CH
  c.hull <- fasterize(hull, pa)
  
  # Mask CH
  c.hull <- mask(c.hull, shp)
  
  dir.create(path = speciesname)
  output_dir <- paste0(speciesname, "/")
  writeRaster(c.hull, paste0(output_dir, "", speciesname, "_CH.tif"), NAflag=-9999)
  
  # Protected area overlap (CH)
  ch.area <- cellStats(c.hull, "sum")
  ov.ch <- cellStats(raster::Which(pa > 0 & c.hull > 0), "sum")
  
  if (!is.finite(ov.ch)) {
    ov.ch <- 0
  }
  
  
  # Saving Model Output
  Output <- 
    data.frame(species = speciesname, 
               ch.area = ch.area, 
               ov.ch = ov.ch)
  
  write.csv(Output, file = paste0("", speciesname, "CH.csv"), row.names = FALSE, quote = TRUE)
  
  # Reduce memory consumption
  dc_cl <- dc_cl[!(dc_cl$species %in% c(i)), ]
  rm("dc", "coords", "hull", "c.hull", "ch.area", "ov.ch", "Output")
  
  
}
