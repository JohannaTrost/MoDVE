library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)  # for arranging plots vertically


indir <- "/Users/johanna/Uni/masterarbeit/data/mc_output/pirineus/scenarios/climdata_era5_cmip6_1906-2024_ssp245_119ts/mc_matrices"

first <- 1906
last <- 2024

# Find mx height across all years
veg_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/pirineus"
max_hgt <- -Inf
for (i in 80:119) {
  veg <- readRDS(file.path(veg_path, paste0("vegp_mof3d_ptm_", i, "_v4.RDS")))
  h <- terra::unwrap(veg$h)
  max_hgt_year <- global(h, "max", na.rm = TRUE)[[1]] + 1
  if (max_hgt_year > max_hgt) {
    max_hgt <- max_hgt_year
  }
}

mc_matrix_ts <- array(NA, dim = c(50, 50, max_hgt, 14, last-first+1))

for (year in first:last) {
  print(year)
  mc_matrix_year <- readRDS(file.path(indir, paste0(year, "_pirineus_mc_matrix.rds")))
  mc_matrix_ts[,, 1:min(dim(mc_matrix_year)[3], max_hgt), , (year-first+1)] <- mc_matrix_year[,, 1:min(dim(mc_matrix_year)[3], max_hgt),]
}

temp_avg <- apply(mc_matrix_ts[,,,1,], 4, mean, na.rm=TRUE)
temp_sd <- apply(mc_matrix_ts[,,,1,], 4, sd, na.rm=TRUE)
rh_avg <- apply(mc_matrix_ts[,,,7,], 4, mean, na.rm=TRUE)
rh_sd <- apply(mc_matrix_ts[,,,7,], 4, sd, na.rm=TRUE)

# --- Plot time series

# Build time vector (1906 to 1950, 45 years)
years <- first:(first + length(first:last) - 1)

# Prepare data frames
df_temp <- data.frame(
  year = years,
  mean = temp_avg,
  sd = temp_sd
)

df_rh <- data.frame(
  year = years,
  mean = rh_avg,
  sd = rh_sd
)

# Plot Temperature
p_temp <- ggplot(df_temp, aes(x = year, y = mean)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              fill = "red", alpha = 0.3) +
  labs(title = "Temperature (1906–2024)",
       y = "Temperature",
       x = "Year") +
  theme_minimal()

# Plot Relative Humidity
p_rh <- ggplot(df_rh, aes(x = year, y = mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              fill = "blue", alpha = 0.3) +
  labs(title = "Relative Humidity (1906–2024)",
       y = "Relative Humidity",
       x = "Year") +
  theme_minimal()

# Arrange vertically
pdf("../../figs/mc_output/temp_1906_2024_pirineus.pdf")
print(p_temp / p_rh)
dev.off()

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
  labs(title = "Temperature by Height (1906–2024)", y = "Temperature", x = "Year") +
  theme_minimal()

p_rh_h <- ggplot(df_rh_h, aes(x = year, y = mean, color = height_bin, fill = height_bin)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(title = "Relative Humidity by Height (1906–2024)", y = "Relative Humidity", x = "Year") +
  theme_minimal()

# Save both plots
pdf("../../figs/mc_output/temp_rh_by_height_1906_2024_pirineus.pdf")
print(p_temp_h / p_rh_h)
dev.off()
