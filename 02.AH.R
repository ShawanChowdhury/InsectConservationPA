options(java.parameters = "-Xmx6g")

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(sf)
library(alphahull)
library(rangeBuilder)
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
  dc <- dc_cl %>% filter(species==i)
  speciesname <- gsub(" ", "_", i)
  dc <- dc %>%
    dplyr::select("decimalLon", "decimalLat")
  
  #Alpha Hull
  range <- getDynamicAlphaHull(dc, coordHeaders=c('decimalLon','decimalLat'), initialAlpha = 2,
                               clipToCoast = 'no', buff = 0, proj = "+proj=longlat +datum=WGS84")
  a.range <- range[[1]]
  alpha_range <- as(a.range, "sf") 
  alpha_range <- st_make_valid(alpha_range)
  alpha_range <- st_transform(alpha_range, crs = wdpa_crs)
  
  # Rasterize AH
  a.hull <- fasterize(alpha_range, pa)
  
  # Mask AH
  a.hull <- mask(a.hull, shp)
  
  dir.create(path = speciesname)
  output_dir <- paste0(speciesname, "/")
  
  writeRaster(a.hull, paste0(output_dir, "", speciesname, "_AH.tif"), NAflag=-9999)
  
  
  # Protected area overlap (AH)
  ah.area <- cellStats(a.hull, "sum")
  ov.ah <- cellStats(raster::Which(pa > 0 & a.hull > 0), "sum")
  
  if (!is.finite(ov.ah)) {
    ov.ah <- 0
  }
  
  
  # Saving Model Output
  Output <- 
    data.frame(species = speciesname, 
               ah.area = ah.area, 
               ov.ah =ov.ah)
  
  write.csv(Output, file = paste0("", speciesname, "AH.csv"), row.names = FALSE, quote = TRUE)
  
}