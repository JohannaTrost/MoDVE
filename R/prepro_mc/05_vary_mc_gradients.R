library(ggplot2)
library(dplyr)

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

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

path <- file.path(base_dir, "climdata_era5_cmip6_1981-2100_ssp245/a1_2")

forest <- "forest0"

mh <- readRDS(file.path(path, forest, "MicrohabitatMatrix199.rds"))

veg_height <- sum(apply(mh[,,, 1], 3, sum) > 0)

mc <- mh[,,, c(3, 4, 5)]
mc[,,, 1] <- mc[,,, 1] * 900
mc_mean <- data.frame(apply(mc, c(3, 4), mean)[1:veg_height,])
names(mc_mean) <- c("mean_light", "mean_relhum", "mean_temp")
mc_sd <- data.frame(apply(mc, c(3, 4), sd)[1:veg_height,])
names(mc_sd) <- c("sd_light", "sd_relhum", "sd_temp")

vertical_mc <- cbind(mc_mean, mc_sd) %>%
  mutate(
    upper_light = mean_light + 1.96 * sd_light,
    upper_relhum = mean_relhum + 1.96 * sd_relhum,
    upper_temp = mean_temp + 1.96 * sd_temp,
    lower_light = mean_light - 1.96 * sd_light,
    lower_relhum = mean_relhum - 1.96 * sd_relhum,
    lower_temp = mean_temp - 1.96 * sd_temp,
  ) %>%
  select(-sd_relhum, -sd_temp) %>%
  pivot_longer(
    cols = -sd_light,
    names_to = c(".value", "Variable"),
    names_pattern = "(mean|upper|lower)_(.*)",
    values_to = c("mean", "upper", "lower")
  ) %>%
  mutate(Variable = case_when(
    Variable == "temp" ~ "Temperature",
    Variable == "relhum" ~ "Relative humidity",
    Variable == "light" ~ "Light",
    TRUE ~ Variable
  )) %>%
  select(-sd_light)
