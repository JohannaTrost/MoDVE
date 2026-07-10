# ------
# Replace vegetation height and PAI in vegp forest data simulated with MoF3D

library(micropoint)
library(microclimf)
library(terra)

# CONFIGURE parameters
REGION <- "regua" # "regua" or "pirineus"
pai_index <- 4 # Index of PAI in microhabitat matrix
min_pai <- 13 # We assume a minimum PAI at every cell with vegetation cover (change if required)
start <- 80 # First time step to consider for vegp from MoF3D
stop <- 199 # Last time step to consider for vegp from MoF3D
ts_plot <- 30 # Plot veg. parameters every 30th time step (set to NA to not generate plots)
rep <- 0 # If no replicate put ""

matrix2raster <- function(matrix, ref_rast, name) {
  out_rast <- copy(ref_rast)
  values(out_rast) <- matrix
  names(out_rast) <- name
  return(out_rast)
}

# Directories
figs_dir <- file.path("..", "modve_figs", "mc_input")
data_dir <- file.path("..", "modve_data")
in_dir <- file.path(data_dir, "mc_input", REGION)
out_dir <- file.path("..", "modve_output")
microhab_dir <- file.path(out_dir, REGION, "microhabitat")

# Verify/create paths
if(!dir.exists(data_dir)) {
  stop(
    paste0(
      data_dir,
      " not found. Make sure vegetation, albedo, and soil data is available (see scripts in 02_surface_inputs)"
    )
  )
}
if(!dir.exists(in_dir)) {
  stop(
    paste0(
      in_dir,
      " not found. Generate forest stand inputs with 02_surface_inputs/04_format_subregion.R."
    )
  )
}
if(!dir.exists(microhab_dir)) {
  stop(
    paste0(
      microhab_dir, " not found. Generate microhabitat matrices from MoF3D (light and forest data) first with",
      " MoDVE/model_pipeline/01_generate_microhabitat.R"
    )
  )
}
if(!dir.exists(figs_dir)) {
  dir.create(figs_dir)
}

# Load veg. parameters derived with the microclimdata package
vegp_reg <- readRDS(paste(in_dir, "vegp.RDS", sep = "/"))

# Unpack 
vegp_unwrpd <- lapply(vegp_reg, terra::unwrap)

# Raster to store simulated forest data
vegp_mof3d <- vegp_unwrpd

# --- Canopy height

for (ts in seq(start, stop)) {

  cat("Processing time step:", ts, "\n")

  # Load MoF3D microhabitat matrix generated with MoDVE/model_pipeline/01_generate_microhabitat.R
  microhab_file <- file.path(microhab_dir, paste0("microhabitatMatrix", ts, ".rds"))
  if(!file.exists(microhab_file)) {
    stop(
      paste0(microhab_file, " not found. Generate microhabitat matrix for time step", ts,
             "with MoDVE/model_pipeline/01_generate_microhabitat.R")
    )
  }
  mm <- readRDS(microhab_file)

  # --- For each cell that is NA extract the max veg. height from PAI

  # Extract PAI
  pai <- mm[,,,pai_index]
  max_x <- dim(pai)[1]
  max_y <- dim(pai)[2]

  # Matrices for canopy height and PAI
  max_heights <- matrix(NA, nrow = max_x, ncol = max_y)
  pai_2d <- matrix(NA, nrow = max_x, ncol = max_y)

  # For each location x,y find max height where PAI > 0
  for (i in 1:max_x) {
    for (j in 1:max_y) {
      # Get the vertical profile of PAI at this cell
      profile <- pai[i, j, ]

      # Find the last layer with non-zero PAI
      non_zero <- which(profile > 0)

      if (length(non_zero) > 0) {
        # Set the height as the maximum non-zero layer (in meters)
        veg_hgt <- max(non_zero)
        max_heights[i, j] <- veg_hgt
        pai_2d[i, j] <- min(min_pai, sum(pai[i, j,]))
      } else {
        max_heights[i, j] <- 0
        pai_2d[i, j] <- 0
      }
    }
  }

  # Store in raster
  vegp_mof3d_ts <- vegp_mof3d
  values(vegp_mof3d_ts$hgt) <- max_heights
  values(vegp_mof3d_ts$pai) <- pai_2d

  if (ts %% ts_plot == 0 & !is.na(ts_plot)) { # Visualize vegetation from model every <ts_plot>-th time step
    print(paste("Plot time step:", ts))

    # Check vegp_mof3d PAI
    pdf(file.path(figs_dir, paste0(REGION, "_veg_mof3d_pai_ts", ts, ".pdf")))
    plot(vegp_mof3d_ts$pai)
    dev.off()

    # Look at average PAI profile
    paii <- apply(pai[,,1:veg_hgt], c(3), mean)
    pdf(file.path(figs_dir, paste0(REGION, "_avg_pai_profile_ts", ts, ".pdf")))
    plot(seq_along(paii) ~ paii, type = "l", main = paste("Total PAI:", sum(paii)))
    dev.off()

    # Check SAI
    sai <- mm[,,,1]
    saii <- apply(sai[,,1:veg_hgt], c(3), mean)
    pdf(file.path(figs_dir, paste0(REGION, "_avg_sai_profile_ts", ts, ".pdf")))
    plot(seq_along(saii) ~ saii, type = "l", main = paste("Total SAI:", sum(saii)))
    dev.off()

    # Check LAI
    lai <- pai - sai
    laii <- apply(lai[,,1:veg_hgt], c(3), mean)
    pdf(file.path(figs_dir, paste0(REGION, "_avg_lai_profile_ts", ts, ".pdf")))
    plot(seq_along(laii) ~ laii, type = "l", main = paste("Total LAI:", sum(laii)))
    dev.off()

    # Check height
    pdf(file.path(figs_dir, paste0(REGION, "_veg_mof3d_hgt_ts", ts, ".pdf")))
    plot(vegp_mof3d_ts$hgt)
    dev.off()
  }

  # Adjust names for micropoint
  names(vegp_mof3d_ts)[names(vegp_mof3d_ts) == "leafr"] <- "lref"
  names(vegp_mof3d_ts)[names(vegp_mof3d_ts) == "leaft"] <- "ltra"
  names(vegp_mof3d_ts)[names(vegp_mof3d_ts) == "hgt"] <- "h"

  # Add coefficient of stomatal conductance sensitivity to photosynthetically active radiation
  vegp_mof3d_ts$q50 <- terra::rast(extent = terra::ext(vegp_mof3d_ts$pai),
                                   resolution = terra::res(vegp_mof3d_ts$pai),
                                   crs = terra::crs(vegp_mof3d_ts$pai))
  vegp_mof3d_ts$q50[] <- 100 # default value if no info is available (see micropoint package "runningmicropoint.Rmd")

  # --- Add veg. emissivity

  em <- 0.97 # From micropoint/microclimf
  vegp_mof3d_ts$em <- deepcopy(vegp_mof3d_ts$q50)
  values(vegp_mof3d_ts$em) <- em

  # Wrap data and save
  vegp_mof3d_wrp <- lapply(vegp_mof3d_ts, terra::wrap)

  # Correct order
  vegp_mof3d_wrp <- vegp_mof3d_wrp[names(micropoint::vegparams)]

  # Save resulting vegp for time step
  saveRDS(vegp_mof3d_wrp, file.path(in_dir, paste0("rep", rep), paste0("vegp_mof3d_ptm_", ts, ".RDS")))
}

