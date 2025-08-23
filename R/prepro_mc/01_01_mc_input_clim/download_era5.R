#!/usr/bin/env Rscript

# ------ Download and process meteorological data - Command Line Version

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript era5_download.R <output_directory> <year>\n")
  cat("Example: Rscript era5_download.R /path/to/output/directory 2042\n")
  quit(status = 1)
}

dirout <- args[1]

# Validate output directory
if (!dir.exists(dirout)) {
  cat("Error: Output directory does not exist:", dirout, "\n")
  cat("Please create the directory first or provide a valid path.\n")
  quit(status = 1)
}

cat("Output directory:", dirout, "\n")

year <- args[2]  # Year to download, e.g., "2024"

# 1983 -> 7-12
# 1984 -> done
# 1985 -> done
# 1986 -> done
# 1987 -> done
# 1988 -> done
# 1989 -> done
# 1990 -> done
# 1993 -> done
# 1994 -> done
# 1995, -> done
# 1997, -> 11-12
# 1998, -> 5-12
# 1999, -> 5-12
# 2000, -> 5-12
# 2001, -> done
# 2002, -> 5-12
# 2003, -> 6-12
# 2004, -> 9-12 missing
# 2005,
# 2007,
# 2009,
# 2010,
# 2011,
# 2012,
# 2013,
# 2014,
# 2015, -> 9-12
# 2018,
# 2019 -> done

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

cat("Loading required libraries...\n")

# Check and install required packages
required_packages <- c("sf", "terra", "mcera5", "ncdf4", "tools", "ggplot2", "patchwork",
                       "foreach", "doParallel", "doRNG", "parallel")

# needs install from github "microclimf", "microclimdata"

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

library(microclimf)
library(microclimdata)

# -- Setup credentials

cat("Setting up credentials...\n")

# Access datasets for download
data(credentials)

uid <- "05d8ef03-8aa5-4508-af48-4c740eeff0d3"
access_token <- "ac0fedd8-e162-4742-949a-a3f40e26a68b"

credentials[3, "username"] <- uid
credentials[3, "password"] <- access_token

# -- Climate Data Setup

cat("Setting up spatial extent and projection...\n")

# Define extent in decimal degrees
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters
r <- project(r, "EPSG:31983", res = 10)

# -- Process year

cat("Downloading", year, "data...\n")

# Time sequence
all_tme <- as.POSIXlt(seq(
      as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
      as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
      by = "1 hour"
))
if (year <= 2003 & year >= 1996) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-05-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 2013) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-04-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 2003) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-07-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 2004) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-09-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 2002) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-11-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 2015) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-12-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else if (year == 1997) {
    tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-11-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
    ))
} else {
  tme <- all_tme
    cat("Time range: Jan 1 - Dec 31,", year, "\n")
}

# Define file prefix
file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")

cat("Downloading ERA5 data...\n")

# Download ERA5 data
req <- era5_download(r, tme, credentials, file_prefix, dirout)


cat("\n=== Downloading complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Downloaded", year, " ERA5 climate data.\n")