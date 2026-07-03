#!/usr/bin/env Rscript

# Merge all cells into one microclimate matrix
# Usage: Rscript model_pipeline/02_simulate_microclimate/03_simulation/02_merge_cells.R --config config.toml [--verbose]

# Example config.toml
#
#     # Dimensions for the microclimate matrix (default: 50x50 if not specified)
#     x_dim = 50
#     y_dim = 50
#
#     # Directories
#     mc_dir = "/path/to/<rep>/<region>/<year>/mc_x<x>_y<y>.rds"  # Directory containing microclimate data
#     veg_dir = "/path/to/<rep>/vegp_mof3d_ptm_<ts>.RDS"   # Directory containing vegetation data
#
#     # Simulation metadata
#     region = "regua"  # "regua" or "pirineus"
#     rep = 1           # Replicate forest number one of {0, 1, 2}
#     year = 2023       # Year of the simulation [1981, 2100]
#     ts = 80           # Time step [80, 199]

suppressPackageStartupMessages({
  library(optparse)
  library(RcppTOML)
  library(terra)
})

# ----------------------------
# CLI argument parsing
# ----------------------------

option_list <- list(
  make_option(c("-c", "--config"), type = "character", help = "Path to TOML config file (no default)"),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  quit(status = 1)
}

# ----------------------------
# Load config
# ----------------------------
tryCatch({
  config <- parseTOML(opt$config)
}, error = function(e) {
  immediateMessage("Error parsing TOML configuration file:", opt$config, "\n")
  immediateMessage("Error message:", e$message, "\n")
  quit(status = 1)
})

x_dim <- if (!is.null(config$x_dim)) config$x_dim else 50
y_dim <- if (!is.null(config$y_dim)) config$y_dim else 50
mc_dir    <- config$mc_dir
veg_dir <- config$veg_dir
region <- config$region
rep <- config$rep
year <- config$year
ts   <- config$ts

# ----------------------------
# Define directories
# ----------------------------

mc_in_dir <- if (!is.na(rep)) file.path(mc_dir, region, paste0("rep", rep), year)
                 else file.path(mc_dir, region, year)

# ----------------------------
# Extract maximum vegetation height
# ----------------------------
vegp_path <- if (!is.na(rep)) file.path(veg_dir, paste0("rep", rep), paste0("vegp_mof3d_ptm_", ts, ".RDS"))
                 else file.path(veg_dir, paste0("vegp_mof3d_ptm_", ts, ".RDS"))
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

