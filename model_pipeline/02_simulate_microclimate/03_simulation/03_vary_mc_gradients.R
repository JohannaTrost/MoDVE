# ------
# Generate vertical microclimate gradient steepness scenarios for sensitivity analysis

library(ggplot2)
library(dplyr)

# ----- GLOBAL variables (to configure)

# start and end time step
start_ts <- 80
end_ts <- 198
# Gradient steepness control
grad_weight <- 1.5 # 0: no gradient, 1: full gradient
# "/path/to/<region>/microhabitat/rep<rep>/microhabitatMatrix<ts>.rds" -> output of 03_add_microclimate_dimensions.R
indir <- file.path("../modve_data/modve_output/regua/scenarios/climdata_era5_cmip6_1906-2024_ssp245_119ts/microhabitat_mc/rep0")
scenario <- paste0("climdata_era5_cmip6_1906-2024_ssp245_119ts_", grad_weight, "_mc_grad")
outdir <- file.path("../modve_data/modve_output/regua/scenarios", scenario, "microhabitat_mc", "rep0")

# -----

dir.create(outdir, recursive = TRUE)

# ----- Flags for microhabitat matrix variables
vars <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
          "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
microhabitatVariableFlags <- c(1, 1, 1, 0, 1, 1, 0)

# ----- Generate gradient steepness scenario

for (ts in seq(start_ts, end_ts)) {
    cat("Time step:", ts, "\n")
    mc <- readRDS(file.path(indir, paste0("microhabitatMatrix", ts, ".rds")))

    mc_no_grad <- mc
    mc_varied_grad <- mc

    humIdx <- which(vars[as.logical(microhabitatVariableFlags)] == "HumNicheOpt")
    tempIdx <- which(vars[as.logical(microhabitatVariableFlags)] == "TempNicheOpt")

    for (mcIdx in c(humIdx, tempIdx)) {
      # Compute vertical means (averaging over z dimension)
      mc_mean <- apply(mc[, , , mcIdx], c(1, 2), mean)

      # Expand back to 3D (50x50x80) by replicating along z
      for (i in 1:dim(mc)[3]) {
        mc_no_grad[, , i, mcIdx] <- mc_mean
      }

      # Compute anomaly
      anom <- mc[, , , mcIdx] - mc_no_grad[, , , mcIdx]

      mc_varied_grad[, , , mcIdx] <- mc_no_grad[, , , mcIdx] + grad_weight * anom

    }

    saveRDS(mc_varied_grad, file.path(outdir, paste0("microhabitatMatrix", ts, ".rds")))
}
