#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(RcppTOML)
  library(terra)
})

# ----------------------------
# CLI argument parsing
# ----------------------------
option_list <- list(
  make_option(c("-c", "--config"), type = "character", help = "Path to TOML config file"),
  make_option(c("-y", "--year"), type = "integer", help = "Year to process"),
  make_option(c("-t", "--ts"), type = "integer", help = "Time step to process")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$year) || is.null(opt$ts)) {
  print_help(opt_parser)
  quit(status = 1)
}

# ----------------------------
# Load config
# ----------------------------
config <- RcppTOML::parseTOML(opt$config)

x_dim     <- config$x_dim
y_dim     <- config$y_dim
mc_dir    <- config$mc_dir
veg_dir <- config$veg_dir

year <- opt$year
ts   <- opt$ts

# ----------------------------
# Define directories
# ----------------------------

mc_in_dir <- file.path(mc_dir, year)

# ----------------------------
# Extract maximum vegetation height
# ----------------------------
vegp_path <- file.path(veg_dir, paste0("vegp_mof3d_ptm_", ts, "_v4.RDS"))
vegp_reg  <- readRDS(vegp_path)
max_hgt   <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE) + 1

n_temp_metrics <- 14
mc_matrix <- array(rep(NA, x_dim * y_dim * max_hgt * n_temp_metrics),
                   dim = c(x_dim, y_dim, max_hgt, n_temp_metrics))

# ----------------------------
# Process cells
# ----------------------------
total_time <- 0
successful_cells <- 0

for (x in 1:x_dim) {
  for (y in 1:y_dim) {
    file_path <- file.path(mc_in_dir, paste0("mc_x", x, "_y", y, ".rds"))

    if (file.exists(file_path)) {
      result <- readRDS(file_path)

      if (!is.null(result)) {
        result$x <- x
        result$y <- y
        actual_heights <- nrow(result$data)

        mc_matrix[x, y, 1:actual_heights, ] <- result$data

        numNAs <- sum(is.na(result$data))
        if (numNAs > 0) {
          cat("Warning: Found", numNAs, "NA values in cell (", x, ",", y, ").\n")
        }

        total_time <- total_time + result$processing_time
        successful_cells <- successful_cells + 1

        if (successful_cells %% 100 == 0) {
          avg_time <- total_time / successful_cells
          cat("Processed", successful_cells, "cells. Average time per cell:",
              round(avg_time, 3), "seconds\n")
        }
      }
    } else {
      stop(paste0("File not found: ", file_path))
    }
  }
}

cat("Processing complete! Total successful cells:", successful_cells, "\n")
cat("Total processing time:", round(total_time, 2), "seconds\n")
cat("Average time per cell:", round(total_time / successful_cells, 3), "seconds\n")

# ----------------------------
# Save result
# ----------------------------
out_file <- file.path(mc_dir, paste(year, region, "mc_matrix.rds", sep = "_"))
saveRDS(mc_matrix, out_file)

cat("Saved result to:", out_file, "\n")

# ----------------------------
# Test reload
# ----------------------------
mc_test <- readRDS(out_file)
