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

install.packages("remotes")
remotes::install_github("dklinges9/mcera5")

library(sf)
library(terra)
library(mcera5)

# example climdata 
library(microclimf)
head(climdata)

# Install microclimdata package from a local path if not yet installed 
remotes::install_github("rspatial/luna") # required for microclimdata
remotes::install_github("ilyamaclean/microclimdata")
library(microclimdata)

# Access datasets for download 
data(credentials)
credentials

# -- Climate Data 

# create template raster
e <- ext(0, 300000, 0, 200000)
r <- rast(e)
crs(r) <- "EPSG:27700"
# create POSIXlt time sequence
tme <- as.POSIXlt(c(0:743)*3600, origin="2022-05-01 00:00", tz = "UTC")
res(r) <- 1000

# Define extent in decimal degrees (WGS84)
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"
res(r) <- 0.25  # Set resolution (in degrees)

# Time sequence for 2024
tme <- as.POSIXlt(seq(
  as.POSIXct("2024-01-01 00:00", tz = "UTC"),
  as.POSIXct("2024-12-31 23:00", tz = "UTC"),
  by = "1 hour"
))

# Define output path and filename
pathout <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/era5_2024"
file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")

# Download ERA5 data
req <- era5_download(r, tme, credentials, file_prefix, pathout)

# ---------- TODO steps below ----------

era5climdata <- era5_process(req, pathout, r, tme)
tc <- rast(era5climdata$temp)
mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", 
                            "orange", "red"))(255)
plot(tc[[1]], col = mypal)
saveRDS(era5climdata, 
        "/Users/johanna/Uni/masterarbeit/code/data/mc_input/processed/era5_climdata_2024.RDS")



