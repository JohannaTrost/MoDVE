
# Load matrix
vars <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
          "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
microhabitatVariableFlags <- c(1, 1, 1, 0, 1, 1, 0)

start_ts <- 80
end_ts <- 198

grad_weight <- 1.5 # 0: no gradient, 1: full gradient

indir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios/climdata_era5_cmip6_1906-2024_ssp245_119ts/a1_2")
scenario <- paste0("climdata_era5_cmip6_1906-2024_ssp245_119ts_", grad_weight, "_mc_grad")
outdir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios", scenario, "a1_2")

dir.create(outdir, recursive = TRUE)

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

# Plot gradients

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

path <- file.path(base_dir, "climdata_era5_cmip6_1906-2024_ssp245_119ts_0.5_mc_grad/a1_2") # TODO plot gradients