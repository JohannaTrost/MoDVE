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

era5_months <- list(
  "1984" = 7:12,
  "1985" = 3:12,
  "1986" = 2:11,
  "1987" = 2:10,
  "1988" = 2:10,
  "1989" = c(3, 4), #8, 9, 10),
  "1990" = c(3, 4), #8, 9, 10),
  "1991" = 2:12,
  "1992" = 2:12,
  "1993" = c(3, 4), # 8, 9, 10),
  "1994" = c(3, 4), #8, 9, 10),
  "1995" = c(5), # 11, 12),
  "1996" = 2:12,
  "1997" = c(5, 6, 7, 8), #12),
  "1998" = c(5), #11, 12),
  "1999" = c(5), #11, 12),
  "2000" = c(5), # 11, 12),
  "2001" = 3:9,
  "2002" = 5,
  "2003" = c(5), #11, 12),
  "2004" = 11,
  "2005" = 4:12,
  "2007" = c(3, 4), #8, 9, 10),
  "2009" = c(3, 4),#8, 9, 10),
  "2010" = c(1), #4, 9, 10),
  "2011" = c(3, 4), #8, 9, 10),
  "2012" = c(3, 4), #8, 9, 10),
  "2013" = c(3, 4), #8, 9, 10),
  "2014" = c(3), #8, 9, 10),
  "2015" = 11,
  "2016" = 4
)

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
if (year == 1984) {
    tme <- as.POSIXlt(seq(as.POSIXct("1984-07-01 00:00", tz="UTC"), as.POSIXct("1984-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1985) {
    tme <- as.POSIXlt(seq(as.POSIXct("1985-03-01 00:00", tz="UTC"), as.POSIXct("1985-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1986) {
    tme <- as.POSIXlt(seq(as.POSIXct("1986-02-01 00:00", tz="UTC"), as.POSIXct("1986-11-30 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1987) {
    tme <- as.POSIXlt(seq(as.POSIXct("1987-02-01 00:00", tz="UTC"), as.POSIXct("1987-10-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1988) {
    tme <- as.POSIXlt(seq(as.POSIXct("1988-02-01 00:00", tz="UTC"), as.POSIXct("1988-10-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1989) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1989-03-01 00:00", tz="UTC"), as.POSIXct("1989-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1989-08-01 00:00", tz="UTC"), as.POSIXct("1989-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1990) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1990-03-01 00:00", tz="UTC"), as.POSIXct("1990-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1990-08-01 00:00", tz="UTC"), as.POSIXct("1990-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1991) {
    tme <- as.POSIXlt(seq(as.POSIXct("1991-02-01 00:00", tz="UTC"), as.POSIXct("1991-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1992) {
    tme <- as.POSIXlt(seq(as.POSIXct("1992-02-01 00:00", tz="UTC"), as.POSIXct("1992-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1993) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1993-03-01 00:00", tz="UTC"), as.POSIXct("1993-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1993-08-01 00:00", tz="UTC"), as.POSIXct("1993-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1994) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1994-03-01 00:00", tz="UTC"), as.POSIXct("1994-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1994-08-01 00:00", tz="UTC"), as.POSIXct("1994-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1995) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1995-05-01 00:00", tz="UTC"), as.POSIXct("1995-05-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1995-11-01 00:00", tz="UTC"), as.POSIXct("1995-12-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1996) {
    tme <- as.POSIXlt(seq(as.POSIXct("1996-02-01 00:00", tz="UTC"), as.POSIXct("1996-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1997) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1997-05-01 00:00", tz="UTC"), as.POSIXct("1997-08-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1997-12-01 00:00", tz="UTC"), as.POSIXct("1997-12-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 1998) {
    tme <- as.POSIXlt(seq(as.POSIXct("1998-05-01 00:00", tz="UTC"), as.POSIXct("1998-05-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 1999) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("1999-05-01 00:00", tz="UTC"), as.POSIXct("1999-05-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("1999-11-01 00:00", tz="UTC"), as.POSIXct("1999-12-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2000) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2000-05-01 00:00", tz="UTC"), as.POSIXct("2000-05-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2000-11-01 00:00", tz="UTC"), as.POSIXct("2000-12-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2001) {
    tme <- as.POSIXlt(seq(as.POSIXct("2001-03-01 00:00", tz="UTC"), as.POSIXct("2001-09-30 23:00", tz="UTC"), by="1 hour"))
} else if (year == 2002) {
    tme <- as.POSIXlt(seq(as.POSIXct("2002-05-01 00:00", tz="UTC"), as.POSIXct("2002-05-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 2003) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2003-05-01 00:00", tz="UTC"), as.POSIXct("2003-05-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2003-11-01 00:00", tz="UTC"), as.POSIXct("2003-12-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2004) {
    tme <- as.POSIXlt(seq(as.POSIXct("2004-11-01 00:00", tz="UTC"), as.POSIXct("2004-11-30 23:00", tz="UTC"), by="1 hour"))
} else if (year == 2005) {
    tme <- as.POSIXlt(seq(as.POSIXct("2005-04-01 00:00", tz="UTC"), as.POSIXct("2005-12-31 23:00", tz="UTC"), by="1 hour"))
} else if (year == 2007) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2007-03-01 00:00", tz="UTC"), as.POSIXct("2007-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2007-08-01 00:00", tz="UTC"), as.POSIXct("2007-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2009) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2009-03-01 00:00", tz="UTC"), as.POSIXct("2009-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2009-08-01 00:00", tz="UTC"), as.POSIXct("2009-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2010) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2010-01-01 00:00", tz="UTC"), as.POSIXct("2010-01-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2010-04-01 00:00", tz="UTC"), as.POSIXct("2010-04-30 23:00", tz="UTC"), by="1 hour"),
        #seq(as.POSIXct("2010-09-01 00:00", tz="UTC"), as.POSIXct("2010-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2011) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2011-03-01 00:00", tz="UTC"), as.POSIXct("2011-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2011-08-01 00:00", tz="UTC"), as.POSIXct("2011-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2012) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2012-03-01 00:00", tz="UTC"), as.POSIXct("2012-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2012-08-01 00:00", tz="UTC"), as.POSIXct("2012-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2013) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2013-03-01 00:00", tz="UTC"), as.POSIXct("2013-04-30 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2013-08-01 00:00", tz="UTC"), as.POSIXct("2013-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2014) {
    tme <- as.POSIXlt(c(
        seq(as.POSIXct("2014-03-01 00:00", tz="UTC"), as.POSIXct("2014-03-31 23:00", tz="UTC"), by="1 hour")
        #seq(as.POSIXct("2014-08-01 00:00", tz="UTC"), as.POSIXct("2014-10-31 23:00", tz="UTC"), by="1 hour")
    ))
} else if (year == 2015) {
    tme <- as.POSIXlt(seq(as.POSIXct("2015-11-01 00:00", tz="UTC"), as.POSIXct("2015-11-30 23:00", tz="UTC"), by="1 hour"))
} else if (year == 2016) {
    tme <- as.POSIXlt(seq(as.POSIXct("2016-04-01 00:00", tz="UTC"), as.POSIXct("2016-04-30 23:00", tz="UTC"), by="1 hour"))
} else {
  tme <- all_tme
    cat("Time range: Jan 1 - Dec 31,", year, "\n")
}

# Define file prefix
file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")

cat("Downloading ERA5 data...\n")

# Download ERA5 data
req <- era5_download(r, tme, credentials, file_prefix, dirout)


# Check files month by month
year_str <- as.character(year)  # current year as string
months_to_check <- era5_months[[year_str]]

for (month in months_to_check) {
  month_str <- sprintf("%02d", month)
  nc_file_path <- file.path(dirout, paste0(file_prefix, year_str, "_", month, ".nc"))

  # Check if dates match
  nc_file <- rast(nc_file_path)
  depth_values <- depth(nc_file)
  dates <- as.POSIXct(depth_values, origin = "1970-01-01", tz = "UTC")
  actual_month <- unique(format(dates, "%m"))

  cat("Expected month:", month_str, "\n")
  cat("Actual month:", actual_month, "\n")
  cat(length(dates), "timestamps\n")
  cat(length(dates) / 24, "days\n")
  cat("----------\n")

  if (length(actual_month) == 1 && actual_month != month_str) {
    cat("⚠️ Mismatch:", month_str, "→", actual_month, "\n\n")
  } else if (length(actual_month) > 1) {
    cat("⚠️ Multiple months detected in", basename(nc_file_path), ":", paste(actual_month, collapse = ", "), "\n\n")
  }
}

cat("\n=== Downloading complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Downloaded", year, " ERA5 climate data.\n")