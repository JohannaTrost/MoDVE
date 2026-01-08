#library(devtools)
#devtools::install("/Users/johanna/Uni/masterarbeit/code/micropoint")
library(micropoint)
library(microclimf)
library(terra)
library(xts)


matrix2raster <- function(matrix, ref_rast, name) {
  out_rast <- copy(ref_rast)
  values(out_rast) <- matrix
  names(out_rast) <- name
  return(out_rast)
}

region <- "pirineus"

# Veg. parameters derived with the microclimdata package
in_dir <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_input/", region)
vegp_reg <- readRDS(paste(in_dir, "vegp_v2.RDS", sep = "/"))

# Unpack 
vegp_unwrpd <- lapply(vegp_reg, terra::unwrap)

# Raster to store simulated forest data
vegp_mof3d <- vegp_unwrpd

# --- Canopy height 

start <- 80
stop <- 99

for (ts in seq(start, stop)) {

  cat("Processing time step:", ts, "\n")

  # MoF3D microhabitat matrix (generated with modified version of A1.R from MoDVE)
  microhab_file <- paste0("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/a1_1/microhabitatMatrix", ts, ".rds")
  mm <- readRDS(microhab_file)

  # --- For each cell that is NA extract the max veg. height from PAI

  # Extract PAI
  pai <- mm[,,,4]
  max_x <- dim(pai)[1]
  max_y <- dim(pai)[2]

  # Create a copy of max_heights to modify
  max_heights <- matrix(NA, nrow = max_x, ncol = max_y)
  pai_2d <- matrix(NA, nrow = max_x, ncol = max_y)

  # For each NA, find max height where PAI > 0
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
        pai_2d[i, j] <- min(13, sum(pai[i, j,]))
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

  if (ts %% 30 == 0) {
    print(paste("Plot time step:", ts))

    # Plot 1: vegp_mof3d PAI
    pdf(paste0("../../figs/mc_input/", region, "_veg_mof3d_pai_ts", ts, ".pdf"))
    plot(vegp_mof3d_ts$pai)
    dev.off()

    # Look at average PAI profile
    paii <- apply(pai[,,1:veg_hgt], c(3), mean)
    pdf(paste0("../../figs/mc_input/", region, "_avg_pai_profile_ts", ts, ".pdf"))
    plot(c(1:length(paii)) ~ paii, type = "l", main = paste("Total PAI:", sum(paii)))
    dev.off()

    # Check SAI
    sai <- mm[,,,1]
    saii <- apply(sai[,,1:veg_hgt], c(3), mean)
    pdf(paste0("../../figs/mc_input/", region, "_avg_sai_profile_ts", ts, ".pdf"))
    plot(c(1:length(saii)) ~ saii, type = "l", main = paste("Total SAI:", sum(saii)))
    dev.off()

    # Check LAI
    lai <- pai - sai
    laii <- apply(lai[,,1:veg_hgt], c(3), mean)
    pdf(paste0("../../figs/mc_input/", region, "_avg_lai_profile_ts", ts, ".pdf"))
    plot(c(1:length(laii)) ~ laii, type = "l", main = paste("Total LAI:", sum(laii)))
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
  vegp_mof3d_ts$q50[] <- 100 # default value if no info is available

  # --- Add veg. emissivity

  em <- 0.97 # From original model (see micropoint/microclimf)
  vegp_mof3d_ts$em <- deepcopy(vegp_mof3d_ts$q50)
  values(vegp_mof3d_ts$em) <- em

  # Wrap data and save
  vegp_mof3d_wrp <- lapply(vegp_mof3d_ts, terra::wrap)

  # Correct order
  vegp_mof3d_wrp <- vegp_mof3d_wrp[names(micropoint::vegparams)]

  saveRDS(vegp_mof3d_wrp, paste0(in_dir, "/vegp_mof3d_ptm_", ts, "_v4.RDS"))

  if (ts %% 30 == 0) {
    # Plot height
    pdf(paste0("../../figs/mc_input/", region, "_veg_mof3d_hgt_ts", ts, ".pdf"))
    plot(vegp_mof3d_ts$h)
    dev.off()
  }
}

