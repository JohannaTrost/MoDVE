require(devtools)
# install_github("ilyamaclean/microclimf")
#remotes::install_local("/Users/johanna/Uni/masterarbeit/code/micropoint",
#                       force = TRUE)

library(micropoint)
library(terra)
library(readr)
library(viridis)


extract_params <- function(raster_list, lon, lat, crs = "EPSG:4326") {
  # Create a SpatVector point in WGS84 (decimal degrees)
  point <- vect(cbind(lon, lat), crs = crs)
  
  # Transform the point to the CRS of the rasters 
  target_crs <- crs(unwrap(raster_list[[1]]))
  point_proj <- project(point, target_crs)
  
  # Extract the value at the projected point for each raster
  result <- lapply(raster_list, function(r) {
    r <- unwrap(r)
    extract(r, point_proj)[[2]]
  })

  return(result)
}

coords2indices <- function(lon, lat, raster, coord_crs = "EPSG:4326") {
  # Create a point in WGS84
  point <- terra::vect(cbind(lon, lat), crs = coord_crs)
  point_proj <- terra::project(point, terra::crs(raster)) # match raster crs
  
  # Get the cell number from the transformed point
  cell <- terra::cellFromXY(raster, 
                            t(as.matrix(terra::geom(point_proj)[, c("x", "y")])))
  
  # Convert the cell number to row and column
  rc <- terra::rowColFromCell(raster, cell)
  names(rc) <- c("x", "y")
  rc <- rc[c("x", "y")]
  
  return(rc)
}

indices2coords <- function(x, y, raster, crs_out = "EPSG:4326") {
  # Get the cell number from row and column indices
  cell <- terra::cellFromRowCol(raster, row = x, col = y)
  
  # Get the center coordinates of the cell in the raster's CRS
  coords <- terra::xyFromCell(raster, cell)
  
  # Create point in specified crs
  pts <- terra::vect(coords, type = "points", crs = terra::crs(raster))
  pts_proj <- terra::project(pts, crs_out)
  
  # Extract coordinates
  coords_proj <- terra::geom(pts_proj)[, c("x", "y")]
  names(coords_proj) <- c("x", "y")
  
  return(coords_proj)
}


# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid 


# Load data for one year 
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp_mof3d_ptm_v2.RDS", sep = "/"))
dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
soilc_reg <- readRDS(paste(in_dir, "soilc_ptm.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))
microhab_file <- "/Users/johanna/Uni/masterarbeit/code/output/modev_zach_25_01_07/MicrohabitatMatrix99.rds"
pai <- readRDS(microhab_file)[,,,5]

# Veg heights
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
min_veg_height <- min(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
n_heights <- length(heights)

# Coordinate 
#lat <- -22.426880
#lon <- -42.765096
indices <- c(1, 1)
coords <- indices2coords(indices[[1]], indices[[2]], terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon <- coords[[1]]
lat <- coords[[2]]

# Get parameters for the point
vegparams <- extract_params(vegp_reg, lon, lat)
grndparams <- extract_params(soilc_reg, lon, lat)
indices <- coords2indices(lon, lat, terra::unwrap(vegp_reg$pai))
#paii <- apply(pai, c(3), "mean")
paii <- pai[indices[[1]], indices[[2]], 1:max(vegparams$h, 0.5)]

mout <- micropoint::runpointmodel(climdata_reg, reqhgt = heights[1], vegparams, 
                                  paii, grndparams, lat = lat, long= lon)

vegparams <- vegparams[names(micropoint::forestparams)]
xx <- plotprofile(climdata_reg, hr = 4091, plotout = "tair", vegparams, 
                  paii = paii,grndparams, lat = lat, long = lon)


# Dimensions
nx <- dim(pai)[1]
ny <- dim(pai)[2]
nz <- dim(pai)[3]  # vertical levels

# Initialize the 3D output array
result <- array(NA, dim = c(nx, ny, nz))

# Unwrap PAI layer for coordinate conversion
unwrapped_pai <- terra::unwrap(vegp_reg$pai)

# Correct parameter order
vegp_reg <- vegp_reg[names(micropoint::forestparams)[1:10]]
soilc_reg <- soilc_reg[names(micropoint::groundparams)]
climdata_reg <- climdata_reg[, names(micropoint::climdata)]

start_time <- Sys.time()

# Loop over all grid indices
for (i in 1:nx) {
  for (j in 1:ny) {
    # Get lon/lat for grid cell
    coords <- indices2coords(i, j, unwrapped_pai)[c("x", "y")]
    lon <- coords[[1]]
    lat <- coords[[2]]
    
    # Extract parameters
    vegparams <- extract_params(vegp_reg, lon, lat)
    grndparams <- extract_params(soilc_reg, lon, lat)
    
    # Get correct index for PAI
    indices <- coords2indices(lon, lat, unwrapped_pai)
    
    # Extract PAI profile for this point
    paii <- pai[indices[[1]], indices[[2]],]
    
    print(length(paii))
    
    # Run the profile simulation
    profile <- tryCatch({
      plotprofile(climdata_reg, hr = 4091, plotout = "tair", vegparams, 
                  paii = paii, grndparams, lat = lat, long = lon)
    }, error = function(e) {
      print(paste("plotprofile failed for grid cell (", i, ",", j, ")"))
      rep(NA, nz)
    })
    
    # Ensure the output is the right length
    profile_length <- length(profile$var)
    if (profile_length > nz) {
      result[i, j, ] <- profile$var[1:nz]
    } else {
      result[i, j, 1:profile_length] <- profile$var
    }
  }
}

end_time <- Sys.time()
time_taken <- round(end_time - start_time, 2)
print(paste("Grid cell took", time_taken, "secs"))


