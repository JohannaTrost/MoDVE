#library(devtools)
#devtools::install("/Users/johanna/Uni/masterarbeit/code/micropoint")
library(micropoint)
library(microclimf)
library(terra)


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
microhab_file <- "/Users/johanna/Uni/masterarbeit/code/output/modev_zach_25_01_07/MicrohabitatMatrix99.rds"
mm <- readRDS(microhab_file)

# --- For each cell that is NA extract the max veg. height from PAI

# Extract PAI
pai <- mm[,,,5]
max_x <- dim(pai)[1]
max_y <- dim(pai)[2]

# Create a copy of max_heights to modify
max_heights <- matrix(NA, nrow = max_x, ncol = max_y)

# For each NA, find max height where PAI > 0
for (i in 1:max_x) {
  for (j in 1:max_y) {
    # Get the vertical profile of PAI at this cell
    profile <- pai[i, j, ]
    
    # Find the last layer with non-zero PAI
    non_zero <- which(profile > 0)
    
    if (length(non_zero) > 0) {
      # Set the height as the maximum non-zero layer (in meters)
      max_heights[i, j] <- max(non_zero)
    }
  }
}

# Replace NAs with 0s -> where no trees
max_heights[is.na(max_heights)] <- 0
max_height <- max(max_heights)

# Store in raster
values(vegp_mof3d$hgt) <- max_heights

# --- PAI

# Canopy as upper 4th of forest -> total PAI
lower_hgt <- floor(max_height * (3/4))
upper_hgt <- ceiling(max_height)
values(vegp_mof3d$pai) <- apply(pai[,,lower_hgt:upper_hgt], 
                                c(1, 2), sum)

plot(vegp_mof3d$pai)
plot(vegp_unwrpd$pai)

paii <- apply(pai[,,1:upper_hgt], c(3), mean)
plot(paii, type = "line")

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

em <- 0.97
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
