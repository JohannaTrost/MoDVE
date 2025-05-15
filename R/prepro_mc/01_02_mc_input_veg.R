library(terra)
library(microclimdata)
library(rgee)
library(luna)

ee_Authenticate()

# - Create a template raster used to define the area of download

# Define extent in decimal degrees
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters 
r <- project(r, "EPSG:31983", res = 10)

GoogleDrivefolder<-"GEE_Exports"
pathtopython<-'/Users/johanna/Uni/ws_24_25/capstone/env/env2/bin/python'
projectname<-'ee-jt281-mc'

# --- Land Cover

# Create a template raster used to define the area of download
lcover_download(r, type = "ESA", year = 2024, GoogleDrivefolder, pathtopython, projectname)
# **** Copy data across to a local drive once complete ***
lcover <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")

plot(lcover)

lcover_4326 <- project(lcover, "EPSG:4326")

plot(lcover_4326)

# --- Vegetation Height

res(r) <- 1
vegheight_download(r, GoogleDrivefolder, pathtopython, projectname)

# Visualize
cheight <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/canopy_height_2020_1m.tif")

plot(cheight)
cheight_4326 <- project(cheight, "EPSG:4326")
plot(cheight_4326)


# --- LAI 

# Empirical
pathout<-"/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/modis_lai/"
#dir.create(pathout)
tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2024-01-01", tz = "UTC")  # Data for January 2023
lai_download(r, tme, reso = 500, pathout, credentials, pathtopython)

# Get LAI
lcover <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")
lcover[lcover == 80] <- NA  # set permanent water bodies to NA'
lcover <- lcover * 0
lcover <- aggregate(lcover, 50, fun = 'mean', na.rm = TRUE)
# Perform mosaic
lai <- lai_mosaic(lcover, pathout, reso = 500, msk = TRUE)

plot(lai)
lai_4326 <- project(lai, "EPSG:4326")
lai_cropped <- crop(lai_4326, e)
plot(lai_cropped)

# --- Derive vegetation parameters 

library(microclimf)

# PAI: Plant area index values represent the combined one sided woody and 
#      green vegetation plant area per unit ground area.

habitats <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")
cheight <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/canopy_height_2020_1m.tif")

habitats[habitats == 10] <- 2   # Forest
habitats[habitats == 20] <- 7   # Shrubland
habitats[habitats == 30] <- 10  # Grassland
habitats[habitats == 40] <- 13  # Cropland
habitats[habitats == 50] <- NA  # Build up 
habitats[habitats == 60] <- NA  # Bare / sparse vegetation
habitats[habitats == 80] <- NA  # Water
habitats[habitats == 90] <- 12  # Herbacious wetland

# Remove color table
coltab(habitats) <- NULL

# Save habitats 
writeRaster(habitats, "/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/habitats.tif", 
            filetype = "GTiff", overwrite = TRUE)

# Below not working properly

habitats <- project(habitats, "EPSG:4326")
cheight <- project(cheight, "EPSG:4326")
lai <- project(lai, "EPSG:4326", res=res(habitats))

# Extract longitude and latitude sequences
coords <- as.data.frame(xyFromCell(habitats, 1:ncell(habitats)))

vegp <- vegpfromhab(habitats, hgts = cheight, pai = lai, 
                    coords$y, coords$x, tme)

plot(terra::unwrap(vegp$pai)[[6]])
plot(terra::unwrap(vegp$leafr))

saveRDS(vegp, 
        "/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/vegp.RDS")



