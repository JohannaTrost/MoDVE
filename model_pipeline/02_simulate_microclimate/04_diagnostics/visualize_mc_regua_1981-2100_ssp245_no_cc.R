# ------
# Visualize MC simulations at REGUA plotting time series for different heights

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)  # for arranging plots vertically

# CONFIGURE directories
# /path/to/<year>_<region>_mc_matrix.rds
indir <- "../modve_data_zenodo/mc_output/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245_no_cc/rep1/mc_matrices"
veg_path <- "../modve_data_zenodo/mc_input/regua/rep1"
out_file <- "../modve_figs/mc_output/temp_rh_by_height_1981–2100_no_cc_regua.pdf"

first <- 1981
last <- 2100

# Find mx height across all years
max_hgt <- -Inf
for (i in 80:119) {
  veg <- readRDS(file.path(veg_path, paste0("vegp_mof3d_ptm_", i, ".RDS")))
  h <- terra::unwrap(veg$h)
  max_hgt_year <- global(h, "max", na.rm = TRUE)[[1]] + 1
  if (max_hgt_year > max_hgt) {
    max_hgt <- max_hgt_year
  }
}

mc_matrix_ts <- array(NA, dim = c(50, 50, max_hgt, 14, last-first+1))
actual_max_hgt <- -Inf

for (year in first:last) {
  mc_matrix_year <- readRDS(file.path(indir, paste0(year, "_regua_mc_matrix.rds")))
  hgt <- min(dim(mc_matrix_year)[3], max_hgt)
  mc_matrix_ts[,, 1:hgt, , (year-first+1)] <- mc_matrix_year[,, 1:min(dim(mc_matrix_year)[3], max_hgt),]

  # COunt NAs per year
  orig_nas <- sum(is.na(mc_matrix_year))
  res_nas <- sum(is.na(mc_matrix_ts[,, 1:hgt, , (year-first+1)]))
  if (orig_nas > 0 || res_nas > 0) {
    print(paste("Year:", year, "Original NAs:", orig_nas, "Result NAs:", res_nas))
  }
  if (hgt > actual_max_hgt) {
    actual_max_hgt <- hgt
  }
}

# Check the number of NAs in the matrix
print(paste("Total NAs in temp:", (sum(is.na(mc_matrix_ts[,,,1,])) / prod(dim(mc_matrix_ts[,,,1,]))) * 100, "%"))
print(paste("Total NAs in relhum:", (sum(is.na(mc_matrix_ts[,,,7,])) / prod(dim(mc_matrix_ts[,,,7,]))) * 100, "%"))

# Build time vector (1906 to 1950, 45 years)
years <- first:(first + length(first:last) - 1)

# --- Define height bins (assuming 3rd dimension = height levels)

n_height <- dim(mc_matrix_ts)[3]
height_bins <- cut(1:n_height, breaks = 3, labels = c("Low", "Medium", "High"))

# --- Function to compute averages for each height bin
get_height_bin_means <- function(var_index) {
  res <- lapply(levels(height_bins), function(bin) {
    h_idx <- which(height_bins == bin)
    mean_vals <- apply(mc_matrix_ts[,,h_idx,var_index,], 4, mean, na.rm = TRUE)
    sd_vals   <- apply(mc_matrix_ts[,,h_idx,var_index,], 4, sd, na.rm = TRUE)
    data.frame(
      year = years,
      mean = mean_vals,
      sd   = sd_vals,
      height_bin = bin
    )
  })
  bind_rows(res)
}

# --- Compute for temperature (var_index = 1) and humidity (var_index = 7)
df_temp_h <- get_height_bin_means(1)
df_rh_h   <- get_height_bin_means(7)

# --- Plots
p_temp_h <- ggplot(df_temp_h, aes(x = year, y = mean, color = height_bin, fill = height_bin)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(title = "Temperature by Height (1981–2100)", y = "Temperature", x = "Year") +
  theme_minimal()

p_rh_h <- ggplot(df_rh_h, aes(x = year, y = mean, color = height_bin, fill = height_bin)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(title = "Relative Humidity by Height (1906–2024)", y = "Relative Humidity", x = "Year") +
  theme_minimal()

# Save both plots
pdf(out_file)
print(p_temp_h / p_rh_h)
dev.off()