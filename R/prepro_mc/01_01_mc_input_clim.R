# ------ Download and process meteorological data 

# obs_time – UTC POSIXlt object of observation times for each climate variable, 2017-01-01 00:00:00
# temp – temperatures (deg C)
# relhum - relative humidity (percentage)
# pres - atmospheric pressure (kPa)
# swdown - total downward shortwave radiation received by a horizontal surface (W/m^2)
# difrad - diffuse radiation (W/m^2)
# lwdown - total downward longward radiation (W/m^2)
# windspeed - wind speed at reference height (m/s)
# winddir - wind direction in degrees
# precip - hourly precipitation (mm). 

# -- Install and load libraries

#install.packages("remotes")
#remotes::install_github("dklinges9/mcera5")

library(sf)
library(terra)
library(mcera5)
library(ncdf4)
library(tools)
library(ggplot2)

# example climdata 
library(microclimf)
head(climdata)

# Install microclimdata package from a local path if not yet installed 
#remotes::install_github("rspatial/luna") # required for microclimdata
#remotes::install_local("/Users/johanna/Uni/masterarbeit/code/microclimdata")
library(microclimdata)

# -- Helper functions

rename_nc_subset <- function(file_path, old_name, new_name) {
  
  tmp_path <- tempfile(fileext = ".nc")
  
  # Open the original NetCDF file
  nc_old <- nc_open(file_path)
  
  # Gather dimensions
  dims <- nc_old$dim
  dim_defs <- lapply(dims, function(d) ncdim_def(d$name, d$units, d$vals))
  
  # Gather variables
  vars <- nc_old$var
  var_defs <- list()
  for (v in names(vars)) {
    var <- vars[[v]]
    var_name <- if (v == old_name) new_name else v
    var_defs[[var_name]] <- ncvar_def(
      name = var_name,
      units = var$units,
      dim = var$dim,
      missval = var$missval,
      longname = var$longname,
      prec = var$prec
    )
  }
  
  # Create new NetCDF file with updated variable name
  nc_new <- nc_create(tmp_path, vars = var_defs)
  
  # Copy over data for each variable
  for (v in names(vars)) {
    old_name_v <- v
    new_name_v <- if (v == old_name) new_name else v
    data <- ncvar_get(nc_old, old_name_v)
    ncvar_put(nc_new, new_name_v, data)
  }
  
  # Close files
  nc_close(nc_old)
  nc_close(nc_new)
  
  # Replace the original file with the modified one
  file.copy(tmp_path, nc_path, overwrite = TRUE)
  unlink(tmp_path)  # Clean up
}

# Access datasets for download 
data(credentials)
credentials


uid <- "05d8ef03-8aa5-4508-af48-4c740eeff0d3"
access_token <- "7dc4c6aa-1201-44e5-bf2a-aa8063aaaba1"

credentials[3, "username"] <- uid
credentials[3, "password"] <- access_token

# -- Climate Data 

# Define extent in decimal degrees
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters 
r <- project(r, "EPSG:31983", res = 10)

years <- setdiff(1981:2023, 2016)
for (year in years) {

  # Time sequence
  if (year == "1981") {
      tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-05-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
      ))
  } else {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
      ))
  }


  # Define output path and filename
  pathout <- "/Users/johanna/Uni/masterarbeit/data/mc_input/climate/era5"
  file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)
                             ),
                        "_")

  # Download ERA5 data
  req <- era5_download(r, tme, credentials, file_prefix, pathout)

  # Recover req object for downloaded files to move on
  req <- build_era5_request(
    xmin = ext(r)$xmin,
    xmax = ext(r)$xmax,
    ymin = ext(r)$ymin,
    ymax = ext(r)$ymax,
    start_time = tme[1],
    end_time = tme[length(tme)],
    outfile_name = file_prefix
  )

  # Correct req object and nc file
  for (month in seq(1, 12)) {
    req[[month]]$target <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)),
                                  paste0("__", year, "_"), month, ".nc")

    # Correct subset name in nc files
    nc_path <- paste0(pathout, "/", req[[month]]$target)
    rename_nc_subset(nc_path, "avg_sdlwrf", "msdwlwrf")
    # avg_sdlwrf -> msdwlwrf
  }

  # Process data
  era5climdata <- era5_process(req, paste(pathout, "/", sep = ""), r, tme)

  # Plot data

  roi <- ext(-42.9, -42.2, -22.8, -22.2)

  # Example palettes
  pal_temp <- colorRampPalette(c("white", "red"))(255)
  pal_relhum <- colorRampPalette(c("white", "blue"))(255)
  pal_pres <- terrain.colors(255)
  pal_swdown <- heat.colors(255)
  pal_difrad <- heat.colors(255)
  pal_lwdown <- heat.colors(255)
  pal_windspeed <- colorRampPalette(c("white", "darkgreen"))(255)
  pal_winddir <- colorRampPalette(c("white", "purple"))(255)
  pal_precip <- colorRampPalette(c("white", "blue", "darkblue"))(255)

  # Variable name to palette mapping
  var_pal <- list(
    temp = pal_temp,
    relhum = pal_relhum,
    pres = pal_pres,
    swdown = pal_swdown,
    difrad = pal_difrad,
    lwdown = pal_lwdown,
    windspeed = pal_windspeed,
    winddir = pal_winddir,
    precip = pal_precip
  )

  # Plotting function
  plot_clim_var_roi <- function(r, title, pal, roi) {
    r_unwr <- unwrap(r)
    r_crop <- crop(r_unwr, roi)
    r_mean <- mean(r_crop, na.rm = TRUE)
    r_df <- as.data.frame(r_mean, xy = TRUE, na.rm = TRUE)
    colnames(r_df)[3] <- "value"

    ggplot(r_df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_gradientn(colours = pal, name = title) +
      coord_equal() +
      theme_minimal() +
      labs(title = title) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  }

  # Loop over all variables
  plots <- lapply(names(era5climdata), function(varname) {
    plot_clim_var_roi(era5climdata[[varname]], varname, var_pal[[varname]], roi)
  })

  # Display plots
  wrap_plots(plots, ncol = 3)

  saveRDS(era5climdata,
          paste0("/Users/johanna/Uni/masterarbeit/data/mc_input/climate/era5_climdata_", year, ".RDS"))
}


