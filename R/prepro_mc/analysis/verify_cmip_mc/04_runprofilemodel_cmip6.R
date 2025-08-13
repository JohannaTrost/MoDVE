#remotes::install_local("/Users/johanna/Uni/masterarbeit/code/forks/micropoint",
#                       force = TRUE)

library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(lubridate)
library(furrr)
library(future)
library(dplyr)

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
    terra::extract(r, point_proj)[[2]]
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

region <- "pirineus"

# Filter relevant dates 2024-09-20 to 2024-09-23 (inclusive)

# start_date <- "2024-10-27 00:00:00" # regua
# end_date <- "2024-10-31 23:59:59" # regua
start_date <- "2024-09-20 00:00:00" # pirineus
end_date <- "2024-09-23 23:59:59" # Pirineus


# Define output directory
outdir <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_output/", region, "_2024_test_v6")
dir.create(outdir)

# Save heights of MC simulations
in_dir <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_input/", region)
vegp_path <- paste(in_dir, "vegp_mof3d_ptm_199_v4.RDS", sep = "/")
vegp_reg <- readRDS(vegp_path)
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)

# Load data for one year
soilc_reg <- readRDS(paste(in_dir, "soilc_v2.RDS", sep = "/"))
climdata_path <- paste(in_dir, "cmip6_climdata_2024_v2.csv", sep = "/")
climdata_reg <- read_csv(climdata_path)
microhab_file <- paste0("/Users/johanna/Uni/masterarbeit/data/modve_output/", region, "/a1/MicrohabitatMatrix199.rds")
pai <- readRDS(microhab_file)[,,,4]

# Veg heights
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height)

climdata_reg <- climdata_reg %>%
  filter(obs_time >= ymd_hms(start_date),
         obs_time <= ymd_hms(end_date))

# Run MC model

# Set grid dimensions
nx <- 4
ny <- 4

# Example: run the model once to get structure
x0 <- 22
y0 <- 22
coords0 <- indices2coords(x0, y0, terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon0 <- coords0[[1]]
lat0 <- coords0[[2]]

vegparams0 <- extract_params(vegp_reg, lon0, lat0)
grndparams0 <- extract_params(soilc_reg, lon0, lat0)
max_hgt0 <- max(values(terra::unwrap(vegp_reg$h)))
paii0 <- pai[x0, y0, 1:max_hgt0]

# Run once to get names, time, height dims
mout_example <- micropoint::runpointmodel(climdata_reg, reqhgt = heights,
                                          vegparams0, paii0, grndparams0,
                                          lat = lat0, long = lon0)

# Initialize 4D arrays for each output variable
output_vars <- names(mout_example)
time_dim <- dim(mout_example[[1]])[2]
height_dim <- dim(mout_example[[1]])[1]

mout_fp_combined <- list()
for (var in output_vars) {
  mout_fp_combined[[var]] <- array(NA_real_, dim = c(nx, ny, height_dim, time_dim))
}

# Main loop over grid
start <- 21
for (i in 1:nx) {
  for (j in 1:ny) {

    x <- start + i
    y <- start + j

    print(paste("Processing grid cell:", x, y))

    coords <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
    lon <- coords[[1]]
    lat <- coords[[2]]

    vegparams <- extract_params(vegp_reg, lon, lat)
    grndparams <- extract_params(soilc_reg, lon, lat)

    max_hgt <- max(values(terra::unwrap(vegp_reg$h)))
    paii <- pai[x, y, 1:max_hgt]

    result <- tryCatch({
      micropoint::runpointmodel(climdata_reg, reqhgt = heights, vegparams,
                                paii, grndparams, lat = lat, long = lon)
    }, error = function(e) {
      warning(sprintf("Model failed at (%d, %d): %s", x, y, e$message))
      NULL
    })

    if (!is.null(result)) {
      for (var in output_vars) {
        mout_fp_combined[[var]][i, j, , ] <- result[[var]]
      }
      print(dim(result$tair))
    }
  }
}

# Creat dir if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
saveRDS(mout_fp_combined, file = paste0(outdir, "/sim_with_era5_",
                                        start + 1 , "_", start + nx, "_",
                                        substr(start_date, 1, 10),
                                        "_",
                                        substr(end_date, 1, 10),
                                        ".rds"))