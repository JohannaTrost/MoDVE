require(devtools)
# install_github("ilyamaclean/microclimf")

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
in_dir <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/regua"
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
indices <- c(30, 25)
coords <- indices2coords(indices[[1]], indices[[2]], terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon <- coords[[1]]
lat <- coords[[2]]

# Get parameters for the point
vegparams <- extract_params(vegp_reg, lon, lat)
grndparams <- extract_params(soilc_reg, lon, lat)
#indices <- coord_2_index(lon, lat, terra::unwrap(vegp_reg$pai))
paii <- apply(pai, c(3), "mean")[1:((max(vegparams$h, 0.5) + 26) * 0.7)]
#pai[indices[[1]], indices[[2]], 1:max(vegparams$h, 0.5)]

mout <- micropoint::runpointmodel(climdata_reg, reqhgt = heights[1], vegparams, 
                                  paii, grndparams, lat = lat, long= lon)


vegp_debug <- micropoint::forestparams
vegp_debug$h <- vegparams$h
vegp_debug$pai <- vegparams$pai
vegp_debug$x <- vegparams$x
vegp_debug$lref <- vegparams$lref
vegp_debug$clump <- vegparams$clump
vegp_debug$ltra <- vegparams$ltra
vegp_debug$leafd <- vegparams$leafd
vegp_debug$em <- vegparams$em
vegp_debug$gsmax <- vegparams$gsmax
vegp_debug$q50 <- vegparams$q50
vegp_debug$skew <- vegparams$skew
vegp_debug$spread <- vegparams$spread

# --- Fit parameters to the data

# Objective function
fit_fun <- function(params, observed, vegparams) {
  skew <- params[1]
  spread <- params[2]
  modeled <- PAIgeometry(PAI = sum(observed), n = length(observed),
                         skew = skew, spread = spread)
  sum((observed - modeled)^2)
}

# Initial guesses for skew and spread
init_params <- c(skew = 1, spread = 1)

# Fit using optim
fit <- optim(par = init_params, fn = fit_fun, observed = paii, vegparams = vegp_debug)

# Fitted parameters
fitted_skew <- fit$par[1]
fitted_spread <- fit$par[2]

paii_test <- PAIgeometry(PAI = sum(paii), 
                         n = length(paii), 
                         skew = fitted_skew, spread = fitted_spread)
z <- c(1:length(paii)) / length(paii)
# plant area within each layer
#plot(z ~ paii_test, type = "l", main = paste("Total PAI:", sum(paii_test)))

#plot(c(1:length(paii)) / length(paii) ~ paii, type = "l", main = paste("Total PAI:", sum(paii)))

vegparams <- vegparams[names(vegp_debug)]
xx <- plotprofile(climdata_reg, hr = 4091, plotout = "tair", vegparams, 
                  paii = paii,grndparams, lat = lat, long = lon)


