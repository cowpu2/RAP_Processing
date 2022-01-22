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
##
## CodeMonkey:  Mike Proctor
## ======================================================================

#library(readr)
library(tidyverse)
library(lubridate)
library(magrittr)
library(rprojroot)
library("tidylog", warn.conflicts = FALSE)
library(raster)
library(sf)
library(gdalUtils)
# 2020-12-16 14:33:00 ------------mdp Incorporate woody transitions,
#                                     and biomass
# 2020-12-17 09:36:37 ------------------------------mdp
# 2021-05-03 10:05:25 ------------------------------mdp Revised specs for NOBO paper


## Local stuff  =================
base_path   <- find_rstudio_root_file()
source_path <- file.path(base_path, "source_data//") ## Shape files for extents go here

##  Loading and setting projection for extents for clipping RAP data -----
#extents_shp <- st_read(paste0(source_path, "/Oswalt_RAP_extents.shp"))
#
crs_str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ##  EPSG:4326 - has to be this for RAP

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
RAP_year <- c(1986:2020) ## a range of years

##  Using version 2 RAP data  - change "v2" to "v3" to access version 3 RAP (I think)
##  Vegetation cover ----
RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v2/vegetation-cover-v2-"
RAP_path
raster_path <- file.path(base_path, "source_data/veg_cover//")
raster_type <- "Cover"

##  Biomass  ------
RAP_path    <- "http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-biomass/v2/vegetation-biomass-v2-"
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
##  Writes rasters to appropriate folder
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


##  Plotting/generate pdfs from all tifs in folder

raster.list <- list.files(raster_path, full.names = TRUE, pattern = ".tif$")
dilbert     <- lapply(raster.list, raster)


##  Just show the plots ------
for (i in seq_along(dilbert)){

    plot(dilbert[[i]],
         main = raster_type,
         sub = dilbert[[i]]@file@name)+
      plot(buffer_shape$geometry, add = TRUE)+
      plot(boundary_shape$geometry, add = TRUE)
}

##  Generate pdfs of plots ----
for (i in seq_along(dilbert)){
  pdf(paste0(dilbert[[i]]@file@name,".pdf"))

  plot(dilbert[[i]],
       main = raster_type,
       sub = dilbert[[i]]@file@name)+
    plot(buffer_shape$geometry, add = TRUE)+
    plot(boundary_shape$geometry, add = TRUE)
  dev.off()
}
