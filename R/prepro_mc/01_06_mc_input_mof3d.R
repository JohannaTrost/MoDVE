#library(devtools)
#devtools::install("/Users/johanna/Uni/masterarbeit/code/micropoint")
library(micropoint)
library(microclimf)
library(terra)
library(xts)


matrix2raster <- function(matrix, ref_rast, name) {
  out_rast <- copy(ref_rast)
  values(out_rast) <- matrix
  names(out_rast) <- name
  return(out_rast)
}


get_by_soiltype <- function(soil_raster) {
  # Soil parameters table (already provided)
  soilparams <- micropoint::soilparams
  
  # Exclude specific columns
  exclude_vars <- c("Soil.type")
  target_vars <- setdiff(colnames(soilparams), exclude_vars)
  
  # soil type codes in raster correspond to rows in soilparams:
  soilparams$soil_code <- 1:nrow(soilparams)
  
  # Create a raster for each variable in soilparams
  raster_list <- list()
  for (var in target_vars) {
    # Create a vector to map soil type code to variable value
    map_values <- setNames(soilparams[[var]], soilparams$soil_code)
    
    # Reclassify the raster
    new_raster <- classify(soil_raster, 
                           rcl = cbind(soilparams$soil_code, soilparams[[var]]))
    names(new_raster) <- var # Set layer name
    
    # Add to list
    raster_list[[var]] <- new_raster
  }
  
  return(raster_list)
}


# Veg. parameters derived with the microclimdata package
in_dir <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp.RDS", sep = "/"))

# Unpack 
vegp_unwrpd <- lapply(vegp_reg, terra::unwrap)

# Raster to store simulated forest data
vegp_mof3d <- vegp_unwrpd

# --- Canopy height 

# MoF3D microhabitat matrix (generated with modified version of A1.R from MoDVE)
microhab_file <- "/Users/johanna/Uni/masterarbeit/code/output/modev_zach_25_01_07/MicrohabitatMatrix98.rds"
mm <- readRDS(microhab_file)

# --- For each cell that is NA extract the max veg. height from PAI

# Extract PAI
pai <- mm[,,,5]
max_x <- dim(pai)[1]
max_y <- dim(pai)[2]

# Create a copy of max_heights to modify
max_heights <- matrix(NA, nrow = max_x, ncol = max_y)
pai_2d <- matrix(NA, nrow = max_x, ncol = max_y)

# For each NA, find max height where PAI > 0
for (i in 1:max_x) {
  for (j in 1:max_y) {
    # Get the vertical profile of PAI at this cell
    profile <- pai[i, j, ]
    
    # Find the last layer with non-zero PAI
    non_zero <- which(profile > 0)
    
    if (length(non_zero) > 0) {
      # Set the height as the maximum non-zero layer (in meters)
      veg_hgt <- max(non_zero)
      max_heights[i, j] <- veg_hgt
      lower_hgt <- floor(veg_hgt * (3/4))  # Get "canopy" PAI only
      pai_2d[i, j] <- sum(pai[i, j, lower_hgt:veg_hgt])
    } else {
      max_heights[i, j] <- 0
      pai_2d[i, j] <- 0
    }
  }
}

# Store in raster
values(vegp_mof3d$hgt) <- max_heights
values(vegp_mof3d$pai) <- pai_2d

plot(vegp_mof3d$pai)
plot(vegp_unwrpd$pai)

# Look at average PAI profile 
paii <- apply(pai[,,1:upper_hgt], c(3), mean)
plot(c(1:length(paii)) ~ paii, type = "l", main = paste("Total PAI:", sum(paii)))

# Check SAI
sai <- mm[,,,1]
saii <- apply(sai[,,1:upper_hgt], c(3), mean)
plot(c(1:length(saii)) ~ saii, type = "l", main = paste("Total SAI:", sum(saii)))

# Check LAI
lai <- pai - sai
laii <- apply(lai[,,1:upper_hgt], c(3), mean)
plot(c(1:length(laii)) ~ laii, type = "l", main = paste("Total LAI:", sum(laii)))

# -- Smoothing - 2D conv

# 2D moving average kernel
kernel <- matrix(1, nrow = 3, ncol = 3)
kernel <- kernel / sum(kernel)

# Apply kernel to every slice
smoothed_pai <- array(0, dim = dim(pai))
for (k in 1:dim(pai)[3]) {
  # prepare cimg object (4D: x, y, cc, z)
  slice <- as.cimg(pai[,,k])
  smoothed <- imfilter(slice, kernel, boundary = "replicate")  # replicate borders
  smoothed_pai[,,k] <- as.array(smoothed)
}

# Adjust names for micropoint 
names(vegp_mof3d)[names(vegp_mof3d) == "leafr"] <- "lref"
names(vegp_mof3d)[names(vegp_mof3d) == "leaft"] <- "ltra"
names(vegp_mof3d)[names(vegp_mof3d) == "hgt"] <- "h"

# Add coefficient of stomatal conductance sensitivity to photosynthetically active radiation
vegp_mof3d$q50 <- terra::rast(extent = terra::ext(vegp_mof3d$pai),
                              resolution = terra::res(vegp_mof3d$pai),
                              crs = terra::crs(vegp_mof3d$pai))
vegp_mof3d$q50[] <- 100 # default value if no info is available

# --- Add veg. emissivity 

em <- 0.97 # From original model (see micropoint/microclimf)
vegp_mof3d$em <- deepcopy(vegp_mof3d$q50)
values(vegp_mof3d$em) <- em

# Wrap data and save 
vegp_mof3d_wrp <- lapply(vegp_mof3d, terra::wrap)

saveRDS(vegp_mof3d_wrp, paste(in_dir, "vegp_mof3d_ptm_v2.RDS", sep = "/"))
saveRDS(paii, paste(in_dir, "paii_mof3d_v2.RDS", sep = "/"))

# Plot height
plot(vegp_mof3d$h)
plot(vegp_unwrpd$hgt)


# --- Get soilparameters for point model 

soilc_reg <- readRDS(paste(in_dir, "soilc.RDS", sep = "/"))

names(soilc_reg)[names(soilc_reg) == "groundr"] <- "gref"

# Extract missing variables
ptm_soilc <- lapply(soilc_reg, terra::unwrap)
soilp <- get_by_soiltype(ptm_soilc$soiltype)

# Add missing variables: Psie, Smax, Smin etc.
for (var_name in names(soilp)) {
  if (var_name == "psi_e") {
    ptm_soilc[[var_name]] <- -soilp[[var_name]]
  } else {
    ptm_soilc[[var_name]] <- soilp[[var_name]]
  }
}
# Rename for correct structure for pointmodel
names(ptm_soilc)[names(ptm_soilc) == "psi_e"] <- "Psie"

dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
ptm_soilc$slope<-terra::terrain(dtm_reg,'slope')
ptm_soilc$aspect<-terra::terrain(dtm_reg,'aspect')

em <- 0.97
ptm_soilc$em <- deepcopy(ptm_soilc$slope)
values(ptm_soilc$em) <- em

ptm_soilc_wrp <- lapply(ptm_soilc, terra::wrap)

saveRDS(ptm_soilc_wrp, paste(in_dir, "soilc_ptm.RDS", sep = "/"))
