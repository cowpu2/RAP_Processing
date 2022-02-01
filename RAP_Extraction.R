## =========== Extracts data from RAP based on extents of shp file=============
##
##  https://rangelands.app
##
##  Includes options to download Vegetative cover, Biomass, and
##  Woody transitions. Woody transitions is available at two resolutions - 30m and 240m
##
##  Specific to Oswalt but can be modified by loading a different extents shape file.
##  This file should be a rectangle/square polygon.  Uses maxx, maxy, minx, miny coordinates of extents file
##
##  Writes rasters to various folders depending on data set
##  Plots and writes pdfs to same folders
##  Files in EPSG4326 folders have EPSG:4326 projection in the file.
##  
##
## CodeMonkey:  Mike Proctor
## ======================================================================

library(readr)
library(tidyverse)
library(lubridate)
library(magrittr)
library(rprojroot)
library("tidylog", warn.conflicts = FALSE)
library(raster)
library(sf)
library(gdalUtils)


## Local stuff  =================
base_path   <- find_rstudio_root_file()
source_path <- file.path(base_path, "source_data//") ## Shape files for extents go here
# raster_path gets defined for each raster type individually below

##  Loading and setting projection for extents for clipping RAP data -----

#crs_str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ##  EPSG:4326 - has to be this for RAP
crs_str <- "EPSG:4326"


Oswalt <- st_read(paste0(source_path, "/boundary.shp"))

# Buffer before reprojecting
Oswalt_1500       <- st_buffer(Oswalt, dist = 1500) # 1500m buffer on Oswalt boundary
Oswalt_1500_WGS84 <- st_transform(Oswalt_1500, crs_str)  ##  Reproject vectors
Oswalt_WGS84      <- st_transform(Oswalt, crs_str)  ##  Reproject vectors

# These are used for plotting - change to match your vectors
boundary_shape <- Oswalt_WGS84      # original boundary
buffer_shape   <- Oswalt_1500_WGS84 # buffered boundary


## gdalwarp will use these values in loop below
extents_shp  <- st_bbox(Oswalt_1500_WGS84)  # if you use st_bbox - te = aoi in loop below
#extents_shp   <- extent(Oswalt_1500_WGS84) # if you use extent - te = c(aoi[1], aoi[3], aoi[2], aoi[4])

aoi <- extents_shp

##  Select years to download  -----
#RAP_year <- c(2019)      ## a single year
RAP_year <- c(2017:2020) ## a range of years

##  There are different versions of coverand biomass data  - change "v2" to "v3" to access version 3 RAP 
##  
##  Vegetation cover ----
#RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v2/vegetation-cover-v2-"
RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-"
RAP_path
raster_path <- file.path(base_path, "source_data/veg_cover//")
raster_type <- "Cover"


##  Biomass  ------
#RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-biomass/v2/vegetation-biomass-v2-"
RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-biomass/v3/vegetation-biomass-v3-"
RAP_path
raster_path <- file.path(base_path, "source_data/biomass//")
raster_type <- "Biomass"


##  Woody transitions - there are multiple resolutions available -------
##  Use appropriate setting below

##  240m resolution ------
RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-derivatives/great-plains/woody-transitions/transitions-240m/great-plains-transitions-240-"
RAP_path
raster_type <- "Woody_trans_240"
raster_path <- file.path(base_path, "source_data/woody_trans/240//")
RAP_year <- c(1990, 2000, 2010, 2015:2019)## This is all years available for woody transitions



##  30m resolution ------
RAP_path <- "http://rangeland.ntsg.umt.edu/data/rap/rap-derivatives/great-plains/woody-transitions/transitions-30m/great-plains-transitions-30-"
RAP_path
raster_type <- "Woody_trans_30"
raster_path <- file.path(base_path, "source_data/woody_trans/30//")
RAP_year <- c(1990, 2000, 2010, 2015:2019)## This is all years available for woody transitions



##  Download and clip RAP data for specified years----
##  Writes rasters to appropriate folder -They are EPSG:4326 but I couldn't get gdalwarp to set them to it.

for (i in RAP_year) {

  RAP_url <- paste0("/vsicurl/",RAP_path, i, ".tif")
  print(RAP_url)
  gdalwarp(srcfile = RAP_url,
           dstfile = paste0(raster_path,"RAP.",raster_type,"_", i,".tif"),
           #te = c(aoi[1], aoi[3], aoi[2], aoi[4]), ##  Can't use te = aoi because they aren't in correct order
           te = aoi, ## Have to use st_bbox() above for this to work
           overwrite = TRUE,
           verbose = TRUE )  
}

##  Reload files, set projection, write back out - choked when I tried to use same name even with overwrite=TRUE====
raster.list <- list.files(raster_path, full.names = TRUE, pattern = ".tif$")
df          <- lapply(raster.list, raster)
df          <- lapply(raster.list, brick)


dir.create(paste0(raster_path, "EPSG4326")) # warning if it exists but doesn't stop anything
EPSG_path <- "EPSG4326/"

for (z in seq_along(df)) {
  projection(df[[z]]) <- "EPSG:4326" # Set projection before writing file
  #print(st_crs(df[[z]]))
  #trimmedName <- trimws( paste(df[[z]]@data@names[z]),"r",whitespace = "[\\h\\v]") #@data@names has a trailing space
  trimmedName[[z]] <- str_sub(df[[z]]@data@names[1], 1, -3)
  print(trimmedName[[z]])
  writeRaster(df[[z]],filename = paste0(raster_path, EPSG_path, trimmedName[[z]],"_4326.tif"), overwrite = TRUE)# EPSG code in filename

}


## Just making sure the projection stuck with the files
raster.list <- list.files(paste0(raster_path, EPSG_path),full.names = TRUE, pattern = ".tif$")
checkCRS    <- lapply(raster.list, raster)  
for (g in seq_along(checkCRS)) {
  print(checkCRS[[g]]@crs@projargs)
  
}



##  Plotting/generate pdfs from all tifs in folder

##  Just show the plots ------
for (i in seq_along(checkCRS)) {

  plot(checkCRS[[i]],
       main = raster_type,
       sub = checkCRS[[i]]@data@names) +
    plot(buffer_shape$geometry, add = TRUE) +
    plot(boundary_shape$geometry, add = TRUE)
  
}


##  Generate pdfs of plots ----
for (i in seq_along(checkCRS)) {
  pdf(paste0(checkCRS[[i]]@file@name,".pdf"))

  plot(checkCRS[[i]],
       main = raster_type,
       sub = checkCRS[[i]]@file@name) +
    plot(buffer_shape$geometry, add = TRUE) +
    plot(boundary_shape$geometry, add = TRUE)
  dev.off()
}
