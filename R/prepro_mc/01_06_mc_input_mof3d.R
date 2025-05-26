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
vegp_mof3d <- copy(vegp_unwrpd)

# --- Canopy height 

forest_file <- '/Users/johanna/Uni/masterarbeit/code/output/mof3d_test/Results/trees_replicate_0_time_step_30.txt'
trees <- read.csv(forest_file, sep = "\t", skip = 8)


# Get voxel positions by ceiling the x and y coordinates
trees$X_voxel <- ceiling(trees$x)
trees$Y_voxel <- ceiling(trees$y)

# Determine matrix dimensions (you may adjust based on known grid size)
tmp <- readLines(forest_file)
max_x <- as.numeric(sub("MaxX\t", "", tmp[4]))
max_y <- as.numeric(sub("MaxY\t", "", tmp[5]))
corridor <- as.numeric(sub("Corridor\t", "", tmp[7]))

# Initialize matrix with NA or -Inf to compare heights
max_heights <- matrix(NA, nrow = max_x, ncol = max_y)

# Loop through data and update matrix with max heights
for (i in seq_len(nrow(trees))) {
  x <- ceiling(trees$x[i] - corridor)
  y <- ceiling(trees$y[i] - corridor)
  h <- trees$height[i]
  
  if (is.na(max_heights[x, y]) || h > max_heights[x, y]) {
    max_heights[x, y] <- h
  }
}
# Replace NAs with 0s -> where no trees
max_heights[is.na(max_heights)] <- 0
max_height <- max(max_heights)

# Store in raster
values(vegp_mof3d$hgt) <- max_heights

# --- PAI

# MoF3D microhabitat matrix (generated with modified version of A1.R from MoDVE)
microhab_file <- "../../output/MoDEV_test_v2/MicrohabitatMatrix30.rds"
mm <- readRDS(microhab_file)

# Extract PAI
pai <- mm[,,,5]

# Canopy as upper 4th of forest -> total PAI
lower_hgt <- floor(max_height * (3/4))
upper_hgt <- ceiling(max_height)
values(vegp_mof3d$pai) <- apply(pai[,,lower_hgt:upper_hgt], 
                                c(1, 2), sum)

plot(vegp_mof3d$pai)
plot(vegp_unwrpd$pai)

paii <- apply(pai[,,1:upper_hgt], c(3), mean)
plot(paii, type = "line")


# --- Leaf angle inclination

zdim <- ceiling(max(max_heights, na.rm = TRUE))
mat_weighted_angle_per_cell <- mm[,,1:zdim,4]

# Convert leaf angle to Cambells leaf angle inclination coefficient
# see https://doi.org/10.1016/0168-1923(88)90057-3
compute_x <- function(theta_deg) {
  if (is.na(theta_deg)) {
    return(NA)
  }
  if (theta_deg <= 57.4) { # >1
    x <- ((1 / theta_deg) - 0.0053) / 0.0103
  }else{
    x <- ((1 / theta_deg) - 0.0107) / 0.0066
  }
  return(x)
}

x_matrix <- apply(mat_weighted_angle_per_cell, c(1, 2, 3), 
                  function(x) compute_x(x))
# Leaf angles with an x-value >10 are unrealistic and will be excluded
angles = copy(mat_weighted_angle_per_cell)
angles[x_matrix > 10] = NA
x_matrix[x_matrix > 10] = NA

# Get average leaf angle inclination coefficients
x_matrix_2d <- apply(x_matrix, c(1, 2), function(x) median(x, na.rm = TRUE))
avg_angles <- apply(angles, c(1, 2), function(x) median(x, na.rm = TRUE))

hist(x_matrix_2d)
hist(avg_angles)
hist(values(terra::unwrap(vegp_reg$x)))

# Store in raster
values(vegp_mof3d$x) <- x_matrix_2d

plot(vegp_mof3d$x)
plot(vegp_unwrpd$x)

# Adjust names for micropoint 
names(vegp_mof3d)[names(vegp_mof3d) == "leafr"] <- "lref"
names(vegp_mof3d)[names(vegp_mof3d) == "leaft"] <- "ltra"
names(vegp_mof3d)[names(vegp_mof3d) == "hgt"] <- "h"

# Add coefficient of stomatal conductance sensitivity to photosynthetically active radiation
vegp_mof3d$q50 <- terra::rast(extent = terra::ext(vegp_mof3d$pai),
                              resolution = terra::res(vegp_mof3d$pai),
                              crs = terra::crs(vegp_mof3d$pai))
vegp_mof3d$q50[] <- 100 # default value if no info is available

# Wrap data and save 
vegp_mof3d_wrp <- lapply(vegp_mof3d, terra::wrap)

saveRDS(vegp_mof3d_wrp, paste(in_dir, "vegp_mof3d.RDS", sep = "/"))
saveRDS(paii, paste(in_dir, "paii_mof3d.RDS", sep = "/"))

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
  ptm_soilc[[var_name]] <- soilp[[var_name]]
}

dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
ptm_soilc$slope<-terra::terrain(dtm_reg,'slope')
ptm_soilc$aspect<-terra::terrain(dtm_reg,'aspect')

#soilc_reg$gref <- microclimf:::.rta(soilc_reg$groundr, length(climdata_reg$obs_time))

