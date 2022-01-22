## ================= Get extents for raster in list  ============
##
##  Used when making raster stacks and for extracting data from RAP
##  so that RAP layers match other rasters/shp files for analysis.  Copy file
##  to current project and change file path to source rasters.
##
## CodeMonkey:  Mike Proctor
## ======================================================================

library(tidyverse)   
library(lubridate)   
library(magrittr)    
library(rprojroot)   
library("tidylog", warn.conflicts = FALSE)
library(raster)
library(sf)
library(tictoc)

# 2020-09-14 13:55:02 ------------------------------mdp

## Local stuff  =================
base_path   <- find_rstudio_root_file()                       
source_path <- file.path(base_path, "source_data//")          
dat_path    <- file.path(base_path, "dat_output//")           
plot_path   <- file.path(base_path, "plots//")                
csv_path    <- file.path(base_path, "csv_output//")           
FocalMean_path    <- file.path(base_path, "FocalMeans//")     

#crs_str <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

##  EPSG:4326 - need this to get extents for gdalwarp/extraction - status 1 error otherwise
crs_str_4326 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

##  Load rasters -----
raster.list <- list.files(FocalMean_path, full.names = TRUE, pattern = ".tif$")

FocalMean_rasters    <- lapply(raster.list, raster)  

##  Set projection -------
FocalMean_rasters    <- lapply(FocalMean_rasters,
                               projectRaster, crs = crs_str_4326) 


##  Use this to check/get extents before rasterstack -------
for (i in seq_along(FocalMean_rasters)){
  print(i)
  print(FocalMean_rasters[[i]]@data@names)
  print(FocalMean_rasters[[i]]@ncols)
  print(FocalMean_rasters[[i]]@nrows)
  print(FocalMean_rasters[[i]]@extent)
  plot(FocalMean_rasters[[i]],
       main = FocalMean_rasters[[i]]@data@names)
  
}

for (z in seq_along(FM_stack)){
  print(z)
  print(FM_stack[[z]]@data@names)
  print(FM_stack[[z]]@ncols)
  print(FM_stack[[z]]@nrows)
  print(FM_stack[[z]]@extent)
  plot(FM_stack[[z]],
       main = FM_stack[[z]]@data@names)
  
}

## All layers need to have the same extents before they can go in a raster stack.

FocalMean_rasters[[7]]@extent  ##  gets extents for layer 7 - whatever layer has greatest extents

# Copy from console to source window and use Alt-leftclick-drag to copy just coords.  
# Select coords with Alt-leftclick-drag in 01_data_extraction.R in Rap_Processing and  
# paste.  Alt-leftclick-drag doesn't work in console.  

xmin       : -97.32063 
xmax       : -97.16268 
ymin       : 33.93385 
ymax       : 34.06129 


