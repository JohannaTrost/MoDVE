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


# Veg. parameters derived with the microclimdata package
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp_v2.RDS", sep = "/"))

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
      pai_2d[i, j] <- min(13, sum(pai[i, j,]))
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
paii <- apply(pai[,,1:veg_hgt], c(3), mean)
plot(c(1:length(paii)) ~ paii, type = "l", main = paste("Total PAI:", sum(paii)))

# Check SAI
sai <- mm[,,,1]
saii <- apply(sai[,,1:veg_hgt], c(3), mean)
plot(c(1:length(saii)) ~ saii, type = "l", main = paste("Total SAI:", sum(saii)))

# Check LAI
lai <- pai - sai
laii <- apply(lai[,,1:veg_hgt], c(3), mean)
plot(c(1:length(laii)) ~ laii, type = "l", main = paste("Total LAI:", sum(laii)))

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

# Correct order
vegp_mof3d_wrp <- vegp_mof3d_wrp[names(micropoint::vegparams)]

saveRDS(vegp_mof3d_wrp, paste(in_dir, "vegp_mof3d_ptm_v3.RDS", sep = "/"))

# Plot height
plot(vegp_mof3d$h)
plot(vegp_unwrpd$hgt)
