library(rgdal)
library(rgeos)
library(raster)

tifs <- list.files(path = "GapRasters_Extent/part",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

s <- raster::stack(tifs)

x <- calc(s, sum)

x[is.na(x[])] <- 0

writeRaster(x, "GapSpecies.tif")