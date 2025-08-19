library(terra)
library(microclimdata)
library(microclimf)
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
lcover <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")

plot(lcover)

lcover_4326 <- project(lcover, "EPSG:4326")

plot(lcover_4326)

# --- Vegetation Height

res(r) <- 1
vegheight_download(r, GoogleDrivefolder, pathtopython, projectname)

# Visualize
cheight <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/canopy_height_2020_1m.tif")

plot(cheight)
cheight_4326 <- project(cheight, "EPSG:4326")
plot(cheight_4326)


# --- LAI 

# Empirical
pathout<-"/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/modis_lai/"
#dir.create(pathout)
tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2024-01-01", tz = "UTC")  # Data for January 2023
lai_download(r, tme, reso = 500, pathout, credentials, pathtopython)

# Get LAI
lcover <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")
lcover[lcover == 80] <- NA  # set permanent water bodies to NA'
lcover <- lcover * 0
lcover <- aggregate(lcover, 50, fun = 'mean', na.rm = TRUE)
# Perform mosaic
lai <- lai_mosaic(lcover, pathout, reso = 500, msk = TRUE)

# Save lai 
writeRaster(lai, "/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/lai_2024_500m.tif", 
            filetype = "GTiff", overwrite = TRUE)

plot(lai)
lai_4326 <- project(lai, "EPSG:4326")
lai_cropped <- crop(lai_4326, e)
plot(lai_cropped)

# --- Derive vegetation parameters 

library(microclimf)

# PAI: Plant area index values represent the combined one sided woody and 
#      green vegetation plant area per unit ground area.

hab <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")
cheight <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/canopy_height_2020_1m.tif")

hab[hab == 10] <- 2   # Forest
hab[hab == 20] <- 7   # Shrubland
hab[hab == 30] <- 10  # Grassland
hab[hab == 40] <- 13  # Cropland
hab[hab == 50] <- NA  # Build up 
hab[hab == 60] <- NA  # Bare / sparse vegetation
hab[hab == 80] <- NA  # Water
hab[hab == 90] <- 12  # Herbacious wetland

# Remove color table
coltab(hab) <- NULL

# Save habitats 
writeRaster(hab, "/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/habitats.tif", 
            filetype = "GTiff", overwrite = TRUE)

# Below not working properly

hab <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/habitats.tif")
cheight <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/canopy_height_2020_1m.tif")
lai <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/lai_2024_500m.tif")

hab <- as.factor(hab)
hab <- project(hab, "EPSG:4326", method="near")
cheight <- project(cheight, "EPSG:4326")
lai <- project(lai, "EPSG:4326", res=res(hab))

# Extract longitude and latitude sequences
coords <- as.data.frame(xyFromCell(hab, 1:ncell(hab)))

tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2024-01-01", tz = "UTC") 
vegp <- vegpfromhab(hab, hgts = cheight, pai = NA, 
                    -22.2, -42.2, tme)

# Define output path
out_dir <- "/Users/johanna/Uni/masterarbeit/figs/mc_input"

# Plot 1: Jan PAI
png(file.path(out_dir, "jan_pai.png"))
plot(terra::unwrap(vegp$pai)[[1]], main = "Jan PAI")
dev.off()

# Plot 2: Vegetation height
png(file.path(out_dir, "vegetation_height.png"))
plot(terra::unwrap(vegp$hgt), main = "Vegetation height")
dev.off()

# Plot 3: Leaf angle coefficient
png(file.path(out_dir, "leaf_angle_coefficient.png"))
plot(terra::unwrap(vegp$x), main = "Leaf angle coefficient")
dev.off()

# Plot 4: Maximum stomatal conductance
png(file.path(out_dir, "max_stomatal_conductance.png"))
plot(terra::unwrap(vegp$gsmax), main = "Maximum stomatal conductance")
dev.off()

# Plot 5: Leaf reflectance
png(file.path(out_dir, "leaf_reflectance.png"))
plot(terra::unwrap(vegp$leafr), main = "Leaf reflectance")
dev.off()

# Plot 6: Canopy clumping factor
png(file.path(out_dir, "canopy_clumping_factor.png"))
plot(terra::unwrap(vegp$clump)[[1]], main = "Canopy clumping factor")
dev.off()

# Plot 7: Leaf diameter
png(file.path(out_dir, "leaf_diameter.png"))
plot(terra::unwrap(vegp$leafd), main = "Leaf diameter")
dev.off()

# Plot 8: Leaf transmittance
png(file.path(out_dir, "leaf_transmittance.png"))
plot(terra::unwrap(vegp$leaft), main = "Leaf transmittance")
dev.off()


saveRDS(vegp, 
        "/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/vegp.RDS")



