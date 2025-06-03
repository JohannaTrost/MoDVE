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
min_veg_height <- min(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
n_heights <- length(heights)

# Coordinate 
#lat <- -22.426880
#lon <- -42.765096
indices <- c(25, 25)
coords <- indices2coords(indices[[1]], indices[[2]], terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon <- coords[[1]]
lat <- coords[[2]]

# Get parameters for the point model
vegparams <- extract_params(vegp_reg, lon, lat)
grndparams <- extract_params(soilc_reg, lon, lat)
indices <- coords2indices(lon, lat, terra::unwrap(vegp_reg$pai))
#paii <- apply(pai, c(3), "mean")
paii <- pai[indices[[1]], indices[[2]], 1:max(vegparams$h, 0.5)]
  
# Format microclimf inputs 
mcf_vegp <- lapply(vegp_mcpoint2mcf(vegp_reg), terra::unwrap)
mcf_soilc <- lapply(soilc_mcpoint2mcf(soilc_reg), terra::unwrap)

# Store all/monthly MC
res <- list()
ptm_res <- list()

# Store profile of mid cell (25, 25)
mcp_profile <- data.frame(height = numeric(0), temp = numeric(0), 
                          relhum = numeric(0), windspeed = numeric(0))
mcf_profile <- data.frame(height = numeric(0), temp = numeric(0), 
                          relhum = numeric(0), windspeed = numeric(0))

# Compute microclimate at different heigths
for (i in 1:length(heights)) {
  
  h <- heights[i]
  print(h)
  
  start_time <- Sys.time()
  
  # Run the micropoint point model to get vertical MC profile
  mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams, 
                                    paii, grndparams, lat = lat, long= lon)
  # Run the microclimf point model to get dfo variables for the grid model
  mcf_ptmout <- microclimf::runpointmodel(climdata_reg, reqhgt = h, dtm_reg, 
                                          mcf_vegp, mcf_soilc)
  
  end_time <- Sys.time()
  print(round(end_time - start_time, 2))
  
  # Format model output to use grid model
  gridm_input <- format_gridm_input(mout, mcf_ptmout$dfo, mcf_ptmout$subs,
                                    climdata_reg, mcf_vegp, lat, lon)
  
  # Save ptm result 
  ptm_res[[as.character(h)]] <- gridm_input$weather
  
  # Replace rel. humidity with microclimf values as micropoint is implausible (extremely dry)
  gridm_input$weather$relhum <- mcf_ptmout$weather$relhum
  
  # Subset the micropoint model output to get hottest and coldest days in each month
  micropoint_mx <- subsetpointmodel(gridm_input, tstep = "month", what = "tmax")
  micropoint_mn <- subsetpointmodel(gridm_input, tstep = "month", what = "tmin")
  
  mout_mx <- runmicro(micropoint_mx, reqhgt = h, mcf_vegp, mcf_soilc, dtm_reg,
                      method = "Cpp")
  mout_mn <- runmicro(micropoint_mn, reqhgt = h, 
                      mcf_vegp, mcf_soilc, dtm_reg,
                      method = "Cpp")
  
  # --- Get monthly min, max and mean
  
  # Time vectors
  obs_time_mn <- micropoint_mn$weather$obs_time
  obs_time_mx <- micropoint_mx$weather$obs_time
  
  # Group time steps by month index
  month_labels_mn <- format(floor_date(obs_time_mn, "month"), "%Y-%m")
  month_labels_mx <- format(floor_date(obs_time_mx, "month"), "%Y-%m")
  if(!all(month_labels_mn == month_labels_mx)) {
    stop("tmin and tmax subsets have an unequal amount of days per month.")
  }
  
  for (var in names(mout_mn)) {
    monthly_mc <- get_monthly_mc_stats(mout_mn[[var]], mout_mx[[var]],
                                       month_labels_mn, month_labels_mx)
    # Populate results list
    if (i > 1) {
      res[[var]][,,i,,] <- monthly_mc
    } else {
      dims <- dim(monthly_mc)
      res[[var]] <- array(NA, dim = c(dims[1], dims[2], n_heights, dims[3], dims[4]))
      res[[var]][,,i,,] <- monthly_mc
    }
  }
  
  # Print avg yearly vars 
  mcf_temp_avg <- mean(res$Tz[indices[[1]], indices[[2]], i,, 3])
  mcf_hum_avg <- mean(res$relhum[indices[[1]], indices[[2]], i,, 3])
  mcf_wind_avg <- mean(res$windspeed[indices[[1]], indices[[2]], i,, 3])
  
  mcp_temp_avg <- mean(ptm_res[[as.character(h)]]$temp)
  mcp_hum_avg <- mean(ptm_res[[as.character(h)]]$relhum)
  mcp_wind_avg <- mean(ptm_res[[as.character(h)]]$windspeed)
  
  # Append to data frames
  mcf_profile <- rbind(mcf_profile, 
                       data.frame(height = h, temp = mcf_temp_avg,
                                  relhum = mcf_hum_avg, windspeed = mcf_wind_avg))
  mcp_profile <- rbind(mcp_profile, 
                       data.frame(height = h, temp = mcp_temp_avg, 
                                  relhum = mcp_hum_avg, windspeed = mcp_wind_avg))
  
}

plot(mcf_profile$height ~ mcf_profile$temp, type="line")
plot(mcp_profile$height ~ mcp_profile$temp, type="line")


