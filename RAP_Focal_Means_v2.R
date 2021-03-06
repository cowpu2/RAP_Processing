## === Focal/Moving window spatial analysis of RAP data for NOBO project  =====
##
## CodeMonkey:  Mike Proctor
##
## RAP data from here https://rangelands.app/
## Previously downloaded by RAP_Extraction.R
##
## Source tif files need to have 4 digit year embedded in filename
##
## https://www.neonscience.org/dc-multiband-rasters-r  is a pretty useful resource
##
## ======================================================================

library(tidyverse)   ##  mutate, transform etc.
library(lubridate)   ##  requried for ts to work
library(magrittr)    ##  %>% - pipe
library(rprojroot)   ##  find_rstudio_root_file()
library("tidylog", warn.conflicts = FALSE)
library(raster)
library(sf)
library(stringr)
library(tictoc)


##  Local stuff  =================
base_path   <- find_rstudio_root_file()
source_path <- file.path(base_path, "source_data//")
cover_path   <- file.path(base_path, "source_data/veg_cover//")
biomass_path   <- file.path(base_path, "source_data/biomass//")
woody30_path   <- file.path(base_path, "source_data/woody_trans/30//")
woody240_path   <- file.path(base_path, "source_data/oody_trans/240//")
dat_path    <- file.path(base_path, "dat_output//")
focal_path    <- file.path(base_path, "source_data/focal_means//")
plot_path   <- file.path(base_path, "plots//")
csv_path    <- file.path(base_path, "csv_output//")
EPSG_path <- "EPSG4326/"


crs_str <- "EPSG:4326"  


##  Load rasters -----
raster.list <- list.files(paste0(cover_path, EPSG_path), full.names = TRUE, pattern = ".tif$")

AFG    <- lapply(raster.list, raster, band = 1)  #  Annual forbs & grasses
BG     <- lapply(raster.list, raster, band = 2)  #  Bare Ground
Litter <- lapply(raster.list, raster, band = 3)  #  Litter
PFG    <- lapply(raster.list, raster, band = 4)  #  Perennial forbs & grasses
Shrubs <- lapply(raster.list, raster, band = 5)  #  Shrubs
Trees  <- lapply(raster.list, raster, band = 6)  #  Trees


## Convert lists to RasterStacks  -----------
AFG_stack    <- stack(AFG)
BG_stack     <- stack(BG)
Litter_stack <- stack(Litter)
PFG_stack    <- stack(PFG)
Shrubs_stack <- stack(Shrubs)
Trees_stack  <- stack(Trees)
toc()

##  Add shrub and tree layers together in new rasterstack  -----------
Woody_stack <- Shrubs_stack ##  Sets up a file with the proper extents, copies the slots
for (i  in 1:nlayers(Shrubs_stack)) {

 df <- sum(Shrubs_stack[[i]], Trees_stack[[i]])
  Woody_stack[[i]] <- df
  names(Woody_stack[[i]]) <-
  # paste0('Woody_', ## Rename with function and year
  #        str_extract(names(Shrubs_stack[[i]]), "\\d\\d\\d\\d"))
    paste0("RAP.","Cover_", str_extract(names(Shrubs_stack[[i]]), "\\d\\d\\d\\d")) # Getting year out of layer
plot(Woody_stack[[i]],
     main = Woody_stack@layers[[i]]@data@names,
     sub = "Oswalt Ranch Woody")
}

##  Add herbaceous layers together  ---------

Herbaceous_stack <- AFG_stack ##  Sets up a file with the proper extents, copies the slots
for (i  in 1:nlayers(AFG_stack)) {

  df <-sum(AFG_stack[[i]], PFG_stack[[i]])
  Herbaceous_stack[[i]] <- df
  names(Herbaceous_stack[[i]]) <-
  paste0("RAP.", "Cover_", str_extract(names(AFG_stack[[i]]), "\\d\\d\\d\\d"))
  plot(Herbaceous_stack[[i]],
       main = Herbaceous_stack@layers[[i]]@data@names,
       sub = "Oswalt Ranch Herbaceous")
}


##  Just to see if extents/crs is right  -------------
for (i in 1:nlayers(AFG_stack)) {
  lapply(AFG_stack@layers[i],        plot)
  lapply(BG_stack@layers[i],         plot)
  lapply(Litter_stack@layers[i],     plot)
  lapply(PFG_stack@layers[i],        plot)
  lapply(Shrubs_stack@layers[i],     plot)
  lapply(Trees_stack@layers[i],      plot)
  lapply(Woody_stack@layers[i],      plot)
  lapply(Herbaceous_stack@layers[i], plot)
}
toc()

##  This is a better way to do above plotting
# animate(AFG_stack, n=1, sub = "Annual Herbaceous")
# animate(BG_stack, n=1, sub = "Bare Ground")



##  Read in shapes  ----
property    <- st_read(paste0(source_path, "Boundary.shp"))
#NOBO_points <- st_read(paste0(source_path, "OR_SurveyPoints.shp"))  ##  NOBO survey points

rebuff      <- st_buffer(property, 1500)   ##  Create buffer vector @ 1500m around property line

#crs_str <- 4326
property    <- st_transform(property,    crs_str)
rebuff      <- st_transform(rebuff,      crs_str)
#NOBO_points <- st_transform(NOBO_points, crs_str)
#NOBO_points <- as_Spatial(NOBO_points)          ##  This has to be class sp rather than sf to be used for extract


plot(rebuff$geometry) +
  plot(property$geometry, add = TRUE) #+
  #plot(NOBO_points$geometry, add = TRUE)

dev.off()

##  Mask to buffer -----
##  Can't mask and crop in one loop because extents change

##  Sets up a file with the proper extents, copies the slots
AFG_stack_mask        <- AFG_stack
BG_stack_mask         <- BG_stack
Litter_stack_mask     <- Litter_stack
PFG_stack_mask        <- PFG_stack
Shrubs_stack_mask     <- Shrubs_stack
Trees_stack_mask      <- Trees_stack
Woody_stack_mask      <- Woody_stack
Herbaceous_stack_mask <- Herbaceous_stack

##  Mask to buffer and plot - plotting is just for drama
tic()
for (x in 1:nlayers(AFG_stack)) {
  df <- mask(AFG_stack[[x]],    rebuff)
    AFG_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(BG_stack[[x]],     rebuff)
    BG_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(Litter_stack[[x]], rebuff)
    Litter_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(PFG_stack[[x]],    rebuff)
    PFG_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(Shrubs_stack[[x]], rebuff)
    Shrubs_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(Trees_stack[[x]],  rebuff)
    Trees_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(Woody_stack[[x]],  rebuff)
    Woody_stack_mask[[x]] <- df
     plot(df)
     print(df)
  df <- mask(Herbaceous_stack[[x]],  rebuff)
     Herbaceous_stack_mask[[x]] <- df
     plot(df)
     print(df)

}
toc()
dev.off()


##  Sets up a file with the proper extents, copies the slots,
##  and converts from RasterLayer to RasterStack so structures align
##  when adding layers
AFG_stack_crop        <- crop(AFG_stack_mask[[1]],         rebuff) %>% stack()
BG_stack_crop         <- crop(BG_stack_mask[[1]],          rebuff) %>% stack()
Litter_stack_crop     <- crop(Litter_stack_mask[[1]],      rebuff) %>% stack()
PFG_stack_crop        <- crop(PFG_stack_mask[[1]],         rebuff) %>% stack()
Shrubs_stack_crop     <- crop(Shrubs_stack_mask[[1]],      rebuff) %>% stack()
Trees_stack_crop      <- crop(Trees_stack_mask[[1]],       rebuff) %>% stack()
Woody_stack_crop      <- crop(Woody_stack_mask[[1]],       rebuff) %>% stack()
Herbaceous_stack_crop <- crop(Herbaceous_stack_mask[[1]],  rebuff) %>% stack()


tic()
for (z in 1:nlayers(AFG_stack_mask)) {
  df <- crop(AFG_stack_mask[[z]], rebuff)
   AFG_stack_crop[[z]] <- df
    plot(AFG_stack_crop[[z]],
       main = AFG_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Annual Herbaceous")
    print(df)
  df <- crop(BG_stack_mask[[z]], rebuff)
   BG_stack_crop[[z]] <- df
    plot(BG_stack_crop[[z]],
       main = BG_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Bare Ground")
    print(df)
  df <- crop(Litter_stack_mask[[z]], rebuff)
   Litter_stack_crop[[z]] <- df
    plot(Litter_stack_crop[[z]],
       main = Litter_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Litter")
    print(df)
  df <- crop(PFG_stack_mask[[z]], rebuff)
  PFG_stack_crop[[z]] <- df
    plot(PFG_stack_crop[[z]],
       main = PFG_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Perennial Herbaceous")
    print(df)
  df <- crop(Shrubs_stack_mask[[z]], rebuff)
  Shrubs_stack_crop[[z]] <- df
    plot(Shrubs_stack_crop[[z]],
       main = Shrubs_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Shrubs")
    print(df)
  df <- crop(Trees_stack_mask[[z]], rebuff)
  Trees_stack_crop[[z]] <- df
    plot(Trees_stack_crop[[z]],
       main = Trees_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Trees")
    print(df)
  df <- crop(Woody_stack_mask[[z]], rebuff)
  Woody_stack_crop[[z]] <- df
    plot(Woody_stack_crop[[z]],
       main = Woody_stack_crop@layers[[z]]@data@names,
       sub = "Oswalt Ranch Woody")
    print(df)
  df <- crop(Herbaceous_stack_mask[[z]], rebuff)
  Herbaceous_stack_crop[[z]] <- df
    plot(Herbaceous_stack_crop[[z]],
         main = Herbaceous_stack_crop@layers[[z]]@data@names,
         sub = "Oswalt Ranch Herbaceous")
    print(df)
}
toc()


for (i in 1:nlayers(AFG_stack_crop)) {  ## Rename layers in each stack

  names(AFG_stack_crop[[i]])         <- paste0("AFG_",AFG_stack_crop@layers[[i]]@data@names)
  names(BG_stack_crop[[i]])          <- paste0("BG_",BG_stack_crop@layers[[i]]@data@names)
  names(Litter_stack_crop[[i]])      <- paste0("LIT_",Litter_stack_crop@layers[[i]]@data@names)
  names(PFG_stack_crop[[i]])         <- paste0("PFG_",PFG_stack_crop@layers[[i]]@data@names)
  names(Shrubs_stack_crop[[i]])      <- paste0("SHRUB_",Shrubs_stack_crop@layers[[i]]@data@names)
  names(Trees_stack_crop[[i]])       <- paste0("TREES_",Trees_stack_crop@layers[[i]]@data@names)
  names(Woody_stack_crop[[i]])       <- paste0("WOODY_",Woody_stack_crop@layers[[i]]@data@names)
  names(Herbaceous_stack_crop[[i]])  <- paste0("HERB_",Herbaceous_stack_crop@layers[[i]]@data@names)

}


names(AFG_stack_crop)
names(BG_stack_crop)
names(Litter_stack_crop)
names(PFG_stack_crop)
names(Shrubs_stack_crop)
names(Trees_stack_crop)
names(Woody_stack_crop)
names(Herbaceous_stack_crop)

crop_list <- c(AFG_stack_crop, BG_stack_crop,Litter_stack_crop,
               PFG_stack_crop, Shrubs_stack_crop, Trees_stack_crop,
               Woody_stack_crop, Herbaceous_stack_crop)

##  Stack of all layers for plotting etc. - Not primary source for layer calculations but does get used there
##
All_Layers <- stack(crop_list)  ##  Stack of stacks - creating stack from list caused problems

tic()
for (i in 1:nlayers(All_Layers)) {
  plot(All_Layers[[i]], main = All_Layers[[i]]@data@names, sub = "All Layers") +
    plot(property$geometry, add = TRUE) +
    plot(rebuff$geometry, add = TRUE) #+
    #plot(NOBO_points$geometry, add = TRUE, col = "red", pch = 19)

}
toc()
dev.off()

nlayers(All_Layers)
names(All_Layers)

##  This isn't useful as it rewrites all of the layer names when reloaded
#writeRaster(All_Layers, filename = paste0(plot_path, "Rap_Oswalt_All.tif"))


##  DF of woody layer to check against orig_values  ----------
Add_results <- round(extract(Woody_stack_mask, NOBO_points),2)##  polygon has to be sp rather than sf
Trees <- round(extract(Trees_stack_mask, NOBO_points),2)
Shrubs <- round(extract(Shrubs_stack_mask, NOBO_points),2)

##  Calculate Focal Means for each layer ----------
##  Intended to do this in a separate file but ALL_Layers
##  loses layer names when saved/loaded.

# Multiple resolutions for focal mean filter - RAP data is 30m resolution
# 3X30  = 90
# 9X30  = 270
# 17X30 = 510
# 35X30 = 1050

ScaleList <- c(3, 9, 17, 35) # Set number of cells for focal mean filter

##  Generate plots of focal mean  ====
focal_mean <- All_Layers  ##  Set extents from existing file



# Needed for extracting to points ~line 359-360
#NOBO_points <- as_Spatial(NOBO_points) ##  This has to be class sp rather than sf to be used for extract

tic()
#ScaleList outside of focal mean
for (z in ScaleList) {# Each resolution
    for (x in 1:nlayers(All_Layers)) { # Each layer
    #for (x in 1:10) { # hacked for video
    focal_mean[[x]] <-
      focal(
        All_Layers[[x]],
        w = matrix(1, nrow = z, ncol = z),
        fun = mean,
        na.rm = FALSE)
    names(focal_mean[[x]]) <- paste0(All_Layers@layers[[x]]@data@names, "_",(z*30),"M")
    df <- focal_mean ##  Useful for animations below - turn off lines 369-370
    df_name <- paste0("Focal_", (z*30))
    assign(df_name,df) 

    projection(df[[x]]) <- "EPSG:4326" # Set projection before writing file
    print(st_crs(df[[x]]))
  
    writeRaster(df[[x]],paste0(focal_path, paste(df@layers[[x]]@data@names,"_Focal.tif")))

    plot(focal_mean[[x]],
         main = focal_mean@layers[[x]]@data@names,
         sub = "Focal Mean")
pdf(paste0(focal_path, focal_mean@layers[[x]]@data@names, ".pdf")) ## Before the plot
    plot(focal_mean[[x]],
         main = focal_mean@layers[[x]]@data@names,
         sub = "Focal Mean")
dev.off()

    }
# Use when extracting values from points - not used when generating rasters  
  # df <- as.data.frame(round(extract(focal_mean, NOBO_points),2))
  # assign(paste0("Focal_", (z*30)),df)

}
toc()

## Just making sure the projection stuck with the files
raster.list <- list.files(focal_path, full.names = TRUE, pattern = ".tif$")
#focalCRS    <- lapply(raster.list, raster, band = 1)  
checkCRS    <- lapply(raster.list, raster)  
for (g in seq_along(checkCRS)) {
  print(checkCRS[[g]]@crs@projargs) # should see "+proj=longlat +datum=WGS84 +no_defs"
  
}


#  More extracting stuff ----------
# rmlist <- c("Focal_1050", "Focal_510", "Focal_90")
# rm(list = (rmlist))

# SiteLabels <- c(1:16) # For rows labels
# 
# Focal_df <- bind_cols(Focal_90, Focal_510, Focal_1050) %>% # Binding side by side
#     mutate(Site = c(1:16))  %>% ## Add site labels
#     dplyr::select(Site, everything())  ## raster package has a select as well

#write_csv(Focal_df, paste0(dat_path, "NOBO_Focal_Means_df", ".dat"))  ##  Write data files ----
#write_csv(Focal_df, paste0(csv_path, "NOBO_Focal_Means_df", ".csv"))

##  Plotting - just for grins and giggles from here down -----
Refocal <-
  Focal_df %>% pivot_longer(-Site, names_to = "Year", values_to = "Avg")

Refocoused <-
  Refocal %>% separate(Year, c("Band","Cover", "Year", "Scale"), sep = "_")


Refocoused$Site <- as.factor(Refocoused$Site)


##  Barplot of focal mean across years ----
ggplot(Refocoused, aes(Year, Avg, col = Year)) +
  geom_boxplot() +
  #geom_point() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(rows = vars(Scale), cols = vars(Band)) +
  ggtitle(paste0("Oswalt Ranch RAP Data - 1986-2020\nFocal Means - All Bands - All Scales"))
  ggsave(paste0(plot_path, "OR_Focal_All.pdf"), width = 14, height = 8, units = "in")

