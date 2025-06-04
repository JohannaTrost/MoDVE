require(devtools)
# install_github("ilyamaclean/microclimf")
#remotes::install_local("/Users/johanna/Uni/masterarbeit/code/micropoint",
#                       force = TRUE)

library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(lubridate)
library(furrr)
library(future)


get_monthly_mc_stats <- function(min_arr, max_arr, month_labs_mn, 
                                 month_labs_mx){
  months <- unique(month_labs_mn)
  # Initialize 4D output array: 50 x 50 x 12 x 3
  result_arr <- array(NA, 
                      dim = c(dim(min_arr)[1:2], length(months), 3))
  
  for (i in seq_along(months)) {
    m <- months[i]
    
    # Get indices for this month
    idx_mn <- which(month_labs_mn == m)
    idx_mx <- which(month_labs_mx == m)
    
    # Compute means for this month
    min_mean <- apply(min_arr[,,idx_mn, drop = FALSE], c(1,2), mean, na.rm = TRUE)
    max_mean <- apply(max_arr[,,idx_mx, drop = FALSE], c(1,2), mean, na.rm = TRUE)
    avg_mean <- (min_mean + max_mean) / 2
    
    result_arr[,,i,1] <- min_mean
    result_arr[,,i,2] <- max_mean
    result_arr[,,i,3] <- avg_mean
  }
  
  return(result_arr)
}


format_gridm_input <- function(mc, dfo, subs, climd, vegd_mcf, lat, lon) {
  in_gridm <- list()
  
  # Extract MC output from micropoint
  weather <- tibble(
    temp      = mc$tair,
    relhum    = mc$relhum,
    windspeed = mc$windspeed,
    obs_time  = mc$obs_time
  )
  
  # Generate weather tibble for grid model
  in_gridm$weather <- climd %>% 
    select(obs_time, pres, swdown, difrad, lwdown, winddir, precip) %>% 
    inner_join(weather, ., by = "obs_time")
  
  
  # Add other variables required by grid model
  in_gridm$dfo <- dfo # Take from microclimf point model 
  in_gridm$Tbz <- NA  # Below soil temperature (not applicable here)
  in_gridm$lat <- lat
  in_gridm$long <- lon
  in_gridm$zref <- max(values(vegd_mcf$hgt))
  in_gridm$subs <- subs
  in_gridm$tmeorig <- mc$obs_time
  in_gridm$matemp <- mean(mc$tair, na.rm = TRUE) # only used for refining below ground temp
  
  class(in_gridm) <- "micropoint"
  
  return(in_gridm)
}


vegp_mcpoint2mcf <- function(veg) {
  names(veg)[names(veg) == "lref"] <- "leafr"
  names(veg)[names(veg) == "ltra"] <- "leaft"
  names(veg)[names(veg) == "h"] <- "hgt"
  
  return(veg[names(microclimf::vegp)])
}


soilc_mcpoint2mcf <- function(soil) {
  names(soil)[names(soil) == "gref"] <- "groundr"
  return(soil[names(microclimf::soilc)])
}


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
# ----------------------------------------------------------
# Define output directory 
outdir <- "/Users/johanna/Uni/masterarbeit/data/mc_output"


# ----------------------------------------------------------
# Set up parallel processing
# ----------------------------------------------------------

# Setup parallel plan (adjust based on your machine)
plan(multisession, workers = parallel::detectCores() - 1)  # or use "multicore" on Linux/macOS

# Create coordinate grid
grid <- expand.grid(x = 3:50, y = 1:50)

# Function to process a single (x, y) cell
process_cell <- function(x, y) {

  library(microclimf)

  # Load data for one year
  in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
  vegp_reg <- readRDS(paste(in_dir, "vegp_mof3d_ptm_v2.RDS", sep = "/"))
  dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
  soilc_reg <- readRDS(paste(in_dir, "soilc.RDS", sep = "/"))
  climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))
  microhab_file <- "/Users/johanna/Uni/masterarbeit/code/output/modev_zach_25_01_07/MicrohabitatMatrix99.rds"
  pai <- readRDS(microhab_file)[,,,5]

  # Veg heights
  max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
  heights <- seq(0.5, max_veg_height + 1)

  message(paste("Processing cell:", x, y))
  start_time_cell <- Sys.time()

  coords <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
  lon <- coords[[1]]
  lat <- coords[[2]]

  vegparams <- extract_params(vegp_reg, lon, lat)
  grndparams <- extract_params(soilc_reg, lon, lat)
  paii <- pai[x, y, 1:max(vegparams$h, 0.5)]

  mcf_vegp <- lapply(vegp_mcpoint2mcf(vegp_reg), terra::unwrap)
  mcf_soilc <- lapply(soilc_mcpoint2mcf(soilc_reg), terra::unwrap)

  res <- list()

  for (i in seq_along(heights)) {

    h <- heights[i]

    message(paste(x, y, "height:", h))
    print(h)

    mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams,
                                      paii, grndparams, lat = lat, long = lon)
    message("Point model run completed.")
    mcf_ptmout <- microclimf::runpointmodel(climdata_reg, reqhgt = h, dtm_reg,
                                            mcf_vegp, mcf_soilc)
    mout$relhum <- mcf_ptmout$weather$relhum

    for (var in c("tair", "tcanopy", "relhum", "windspeed")) {
      message(paste(x, y, var, h))
      if (i == 1) {
        res[[var]] <- array(NA, dim = c(length(heights), length(mout[[var]])))
      }
      res[[var]][i, ] <- mout[[var]]
    }
  }

  saveRDS(res, paste0(outdir, "/mc_x", x, "_y", y, "_v1.rds"))
  end_time_cell <- Sys.time()
  message(paste("Finished cell:", x, y, "in", round(end_time_cell - start_time_cell, 2), "seconds"))

  TRUE
}

# Run in parallel
future_map2(grid$x, grid$y, process_cell, .progress = TRUE, .options = furrr_options(seed = 420))

Sys.time()