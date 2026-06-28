# Generate LAI and vegetation input data based on habitat types
# - Based on Microclimf and microclimdata repositories
# src:
#     - https://github.com/ilyamaclean/microclimdata/blob/main/vignettes/instructions.Rmd
#     - https://github.com/ilyamaclean/microclimf/blob/main/vignettes/running-microclimf.Rmd

library(terra)
library(microclimdata)
library(microclimf)
library(rgee)
library(luna)

ee_Authenticate()

# Parameters to configure for using Google Earth Engine (for details see sources above)
GoogleDrivefolder <- NA # e.g. "GEE_Exports"
pathtopython<-NA # path to your python environment e.g. '~/bin/python'
projectname<-NA # e.g. 'ee-x-mc'

if(is.na(GoogleDrivefolder) | is.na(pathtopython) | is.na(projectname)) {
  stop(
    "Please set GoogleDrivefolder, pathtopython, projectname for the Google Erath Engine setup
    (see microclimdata R package)."
  )
}

# Directories
figs_dir <- file.path("..", "modve_figs", "mc_input")
data_dir <- file.path("..", "modve_data")

if (!dir.exists(file.path(data_dir, "mc_input", "vegetation"))) {
  dir.create(file.path(data_dir, "mc_input", "vegetation"), recursive = TRUE)
}

# - Create a template raster used to define the area of download

# Define extent in decimal degrees of respective region in Atlantic forest (Brazil)
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters 
r <- project(r, "EPSG:31983", res = 10)

# -------------------- Download required products -------------------- #

# Download landcover
lcover_download(r, type = "ESA", year = 2024, GoogleDrivefolder, pathtopython, projectname)

# Download Vegetation Height
res(r) <- 1
vegheight_download(r, GoogleDrivefolder, pathtopython, projectname)

# Copy files to following directories
lcover_path <- file.path(data_dir, "mc_input", "vegetation", "ESA_WorldCover_2023_10m.tif")
cheight_path <- file.path(data_dir, "mc_input", "vegetation", "canopy_height_2020_1m.tif")

if(!file.exists(lcover_path) | !file.exists(cheight_path)) {
  stop(
    paste0(
      "Copy downloaded land cover tif file to ", lcover_path, " and canopy height tif file to " , cheight_path,
      " (also available on Zenodo)." # TODO add zenodo link
    )
  )
}

# - LAI

# Empirical
pathout <- file.path(data_dir, "mc_input", "vegetation", "modis_lai")

if (!dir.exists(pathout)) {
  dir.create(pathout, recursive = TRUE)
}

# Estimate LAI from landcover time and MODIS (as lai_download from microclimdata not working here)
lcover <- rast(lcover_path)
lcover[lcover == 80] <- NA  # set permanent water bodies to NA'
lcover <- lcover * 0
lcover <- aggregate(lcover, 50, fun = 'mean', na.rm = TRUE)
# Perform mosaic
lai <- lai_mosaic(lcover, pathout, reso = 500, msk = TRUE)

# Save lai
lai_path <- file.path(data_dir, "mc_input", "vegetation", "lai_2024_500m.tif")
writeRaster(lai, lai_path, filetype = "GTiff")

# --- Derive vegetation parameters 

hab <- rast(lcover_path)

# Convert from ESA landcover types to required habitat types
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
hab_path <- file.path(data_dir, "mc_input", "vegetation", "habitats.tif")
writeRaster(hab, hab_path, filetype = "GTiff", overwrite = TRUE)

# Load LAI and canopy height
lai <- rast(lai_path)
cheight <- rast(cheight_path)

# Convert to same resolution and crs
hab <- as.factor(hab)
hab <- project(hab, "EPSG:4326", method="near")
cheight <- project(cheight, "EPSG:4326")
lai <- project(lai, "EPSG:4326", res=res(hab))

# Extract longitude and latitude sequences
coords <- as.data.frame(xyFromCell(hab, 1:ncell(hab)))

# Get vegp objects
tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2024-01-01", tz = "UTC") 
vegp <- vegpfromhab(hab, hgts = cheight, pai = NA, -22.2, -42.2, tme)

# Save vegetation parameters
saveRDS(vegp,file.path(data_dir, "mc_input", "vegetation", "vegp.RDS"))

# -- Plot the resulting vegetation data

# Plot 1: Jan PAI
png(file.path(figs_dir, "01_jan_pai.png"))
plot(terra::unwrap(vegp$pai)[[1]], main = "Jan PAI")
dev.off()

# Plot 2: Vegetation height
png(file.path(figs_dir, "02_vegetation_height.png"))
plot(terra::unwrap(vegp$hgt), main = "Vegetation height")
dev.off()

# Plot 3: Leaf angle coefficient
png(file.path(figs_dir, "03_leaf_angle_coefficient.png"))
plot(terra::unwrap(vegp$x), main = "Leaf angle coefficient")
dev.off()

# Plot 4: Maximum stomatal conductance
png(file.path(figs_dir, "04_max_stomatal_conductance.png"))
plot(terra::unwrap(vegp$gsmax), main = "Maximum stomatal conductance")
dev.off()

# Plot 5: Leaf reflectance
png(file.path(figs_dir, "05_leaf_reflectance.png"))
plot(terra::unwrap(vegp$leafr), main = "Leaf reflectance")
dev.off()

# Plot 6: Canopy clumping factor
png(file.path(figs_dir, "06_canopy_clumping_factor.png"))
plot(terra::unwrap(vegp$clump)[[1]], main = "Canopy clumping factor")
dev.off()

# Plot 7: Leaf diameter
png(file.path(figs_dir, "07_leaf_diameter.png"))
plot(terra::unwrap(vegp$leafd), main = "Leaf diameter")
dev.off()

# Plot 8: Leaf transmittance
png(file.path(figs_dir, "08_leaf_transmittance.png"))
plot(terra::unwrap(vegp$leaft), main = "Leaf transmittance")
dev.off()



