library(rgdal)
library(rgeos)
library(raster)
library(dplyr)

wdpa_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

pa <- raster("PA_1km_up.tif")
crs(pa) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

tifs <- list.files(path = "GapSpecies/Rasters",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "_AH")]


tifs1 <- as.data.frame(tifs)

colnames(tifs1) <- c("species")

species <- unique(tifs1$species)


for (i in species) {
  x <- raster(i)
  crs(x) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
  
  z <- resample(x, pa, method = "ngb")
  
  z[is.na(z[])] <- 0
  
  crs(z) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
  
  
  sp.name <- gsub("GapSpecies/", "", i)
  sp.name <- as.data.frame(sp.name)
  
  sp.name <- sp.name %>%
    tidyr::separate(sp.name, c("a", "b", "c"), sep = "/")
  
  spname <- sp.name$a
  
  writeRaster(z, paste0("", spname, ".tif"), NAflag=-9999, overwrite = TRUE)
}
