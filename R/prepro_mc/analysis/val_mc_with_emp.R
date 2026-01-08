
library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(lubridate)
library(furrr)
library(future)
library(dplyr)
library(ggplot2)
library(purrr)
library(readxl)
library(tidyverse)
library(patchwork)


extract_params <- function(raster_list, lon, lat, crs = "EPSG:4326") {
  # Create a SpatVector point in WGS84 (decimal degrees)
  point <- vect(cbind(lon, lat), crs = crs)

  # Transform the point to the CRS of the rasters
  target_crs <- crs(unwrap(raster_list[[1]]))
  point_proj <- project(point, target_crs)

  # Extract the value at the projected point for each raster
  result <- lapply(raster_list, function(r) {
    r <- unwrap(r)
    terra::extract(r, point_proj)[[2]]
  })

  return(result)
}

indices2coords <- function(x, y, raster, crs_out = "EPSG:4326") {
  # Get the cell number from row and column indices
  cell <- terra::cellFromRowCol(raster, row = x, col = y)

  # Get the center coordinates of the cell in the raster's CRS
  coords <- terra::xyFromCell(raster, cell)

  # Create point in specified crs
  pts <- terra::vect(coords, type = "points", crs = terra::crs(raster))
  pts_proj <- terra::project(pts, crs_out)

  # Extract coordinates
  coords_proj <- terra::geom(pts_proj)[, c("x", "y")]
  names(coords_proj) <- c("x", "y")

  return(coords_proj)
}

# --- Prep data

# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid

# Function to process a single (x, y) cell
x <- 25
y <- 25

# Load data for one year
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
in_dir_regua <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir_regua, "vegp_mof3d_ptm_v2.RDS", sep = "/"))  # PAI and canopy height will be replaced by MoF3D output below
soilc_reg <- readRDS(paste(in_dir, "soilc_v2.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "climdata_era5_cmip6_2024_v3.csv", sep = "/")) # REGUA
#climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/")) # Pirineaus SA

# Get coordiantes
coords_veg <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
coords_soil <- indices2coords(x, y, terra::unwrap(soilc_reg$aspect))[c("x", "y")]
lon <- coords_soil[[1]]
lat <- coords_soil[[2]]

# Get params for the point
vegparams <- extract_params(vegp_reg, coords_veg[[1]], coords_veg[[2]])
grndparams <- extract_params(soilc_reg, lon, lat)

# Get PAI
microhab_file <- "/Users/johanna/Uni/masterarbeit/data/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/a1_2/forest0/microhabitatMatrix123.rds" # Regua
#microhab_file <- "/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios/climdata_era5_cmip6_1906-2024_ssp245_119ts/a1_2/microhabitatMatrix198.rds" # Pirineus

pai <- readRDS(microhab_file)[,,,5]
paii <- apply(pai[,,1:max(vegparams$h, 0.5)], c(3), mean, na.rm = TRUE)
#paii <- pai[25, 25, 1:max(vegparams$h, 0.5)]

# Replace PAI
vegparams$pai <- sum(paii)

# --- Simulate microclimate

# Veg heights
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
i <- 2
h <- heights[i]

mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams,
                                  paii, grndparams, lat = lat, long = lon)

# put into dataframe
mc_sim <- data.frame(
  obs_time = mout$obs_time,
  tair = mout$tair,
  tleaf = mout$tleaf,
  relhum = mout$relhum
)

# --- Load empirical data

# -- Pirineus
# Define paths
# emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 770m elevation low Pirineus understorey transect"
#
# # Define logger directories
# emp_dirs <- c("10m", "30m", "50m", "70m", "90m")
# macro_dir <- "Reference or at 30m in rockfall sort of gap on pole"  # This is the macroclimate reference

# # -- Regua
# Define paths
emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 400m elevation REGUA understory Trilha Verde"

# Define logger directories
emp_dirs <- c("3600m, 446mASL", "3650m, 438mASL", "3700m, 433mASL", "3750m, 398mASL", "3800m, 387mASL")
macro_dir <- "3850m, 387mASL waterfall reference"  # This is the macroclimate reference

# Function to load and process a logger file
load_logger_data <- function(dir_name) {
  file_path <- file.path(emp_path, dir_name, "mc_data.xlsx")

  read_excel(file_path) %>%
    select("Date-Time (Brazil Standard Time)", "Temperature (°C)", "RH (%)") %>%
    rename(obs_time = "Date-Time (Brazil Standard Time)",
           tair = "Temperature (°C)",
           relhum = "RH (%)") %>%
    mutate(
      obs_time = force_tz(as.POSIXct(obs_time), tzone = "America/Sao_Paulo"),
      obs_time_utc = with_tz(obs_time, tzone = "UTC"),
      obs_hour_utc = floor_date(obs_time_utc, unit = "hour")
    ) %>%
    group_by(obs_hour_utc) %>%
    summarise(
      tair = mean(tair, na.rm = TRUE),
      relhum = mean(relhum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(obs_time = obs_hour_utc)
}

# Load macroclimate data
macro_data <- load_logger_data(macro_dir) %>%
  rename(tair_macro = tair, relhum_macro = relhum)

# --- Compare empirical data against model and macroclimate

# Function to compare one logger against model and macroclimate
process_comparison <- function(dir_name) {
  emp_data <- load_logger_data(dir_name) %>%
    rename(tair_emp = tair, relhum_emp = relhum)

  # Join with model and macro
  joined <- emp_data %>%
    left_join(mc_sim, by = "obs_time") %>%
    left_join(macro_data, by = "obs_time")

  tibble(
    logger = dir_name,

    # Model vs Empirical
    rmse_tair_model = sqrt(mean((joined$tair_emp - joined$tair)^2, na.rm = TRUE)),
    cor_tair_model = cor(joined$tair_emp, joined$tair, use = "complete.obs"),
    rmse_relhum_model = sqrt(mean((joined$relhum_emp - joined$relhum)^2, na.rm = TRUE)),
    cor_relhum_model = cor(joined$relhum_emp, joined$relhum, use = "complete.obs"),

    # Macroclimate vs Empirical (skip if same as macro reference)
    rmse_tair_macro = if (dir_name != macro_dir) sqrt(mean((joined$tair_emp - joined$tair_macro)^2, na.rm = TRUE)) else NA_real_,
    cor_tair_macro = if (dir_name != macro_dir) cor(joined$tair_emp, joined$tair_macro, use = "complete.obs") else NA_real_,
    rmse_relhum_macro = if (dir_name != macro_dir) sqrt(mean((joined$relhum_emp - joined$relhum_macro)^2, na.rm = TRUE)) else NA_real_,
    cor_relhum_macro = if (dir_name != macro_dir) cor(joined$relhum_emp, joined$relhum_macro, use = "complete.obs") else NA_real_
  )
}

# Apply to all loggers
results <- map_dfr(emp_dirs, process_comparison)

# View sorted by model RH MAE
results_sorted <- results %>% arrange(rmse_relhum_model)
print(results_sorted)

# --- Show results of the comparison (MAE and correlation)

# Summarize across loggers (excluding macro itself)
summary_stats <- results %>%
  filter(logger != macro_dir) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), na.rm = TRUE)) %>%
  t() %>%
  round(2)

print(summary_stats)

# --- Some plots

# Load and process all logger data
emp_data_list <- map(emp_dirs, load_logger_data)

# Combine all into a single tibble by full_join by time
repl_tair <- paste0("tair_", emp_dirs[length(emp_dirs)])
repl_relhum <- paste0("relhum_", emp_dirs[length(emp_dirs)])

emp_data_combined <- emp_data_list %>%
  reduce(full_join, by = "obs_time") %>%
  rename_at(vars(starts_with("tair")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename_at(vars(starts_with("relhum")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename(!!repl_tair := "tair", !!repl_relhum := "relhum")

names(emp_data_combined)

# Get era5 / cmip6 macroclimate data
model_macro <- climdata_reg %>%
    select(obs_time, temp, relhum) %>%
    rename(tair_macro_sim = temp, relhum_macro_sim = relhum)

# Compute mean and standard deviation across loggers for each time point
emp_summary <- emp_data_combined %>%
  rowwise() %>%
  mutate(
    tair_mean = mean(c_across(starts_with("tair_")), na.rm = TRUE),
    tair_sd = sd(c_across(starts_with("tair_")), na.rm = TRUE),
    relhum_mean = mean(c_across(starts_with("relhum_")), na.rm = TRUE),
    relhum_sd = sd(c_across(starts_with("relhum_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(obs_time, tair_mean, tair_sd, relhum_mean, relhum_sd) %>%
  inner_join(., macro_data, by = "obs_time") %>%
  inner_join(., mc_sim, by = "obs_time") %>%
  inner_join(., model_macro, by = "obs_time")

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = tair_mean - 1.96 * tair_sd,
                  ymax = tair_mean + 1.96 * tair_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = tair_mean, color = "Measured microclimate", linetype = "Measured microclimate"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = tair_macro,
                color = "Measured macroclimate",
                linetype = "Measured macroclimate"),
            size = 0.6) +
  geom_line(aes(x = obs_time, y = tair, color = "Simulated microclimate", linetype = "Simulated microclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = tair_macro_sim,
                color = "CMIP6 macroclimate",
                linetype = "CMIP6 macroclimate"),
            size = 0.6) +
  scale_color_manual(name = "",
                     values = c("Measured microclimate" = "#3F826D",
                                "Measured macroclimate" = "#3F826D",
                                "Simulated microclimate" = "#FAC05E",
                                "CMIP6 macroclimate" = "#FAC05E")) +

  scale_linetype_manual(name = "",
                        values = c("Measured microclimate" = "solid",
                                   "Measured macroclimate" = "dashed",
                                   "Simulated microclimate" = "solid",
                                   "CMIP6 macroclimate" = "dashed")) +
  labs(
       x = "",
       y = "Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        legend.title = element_text(face = "bold"))

# Second plot: Relative Humidity (empirical vs. simulated)
plot_relhum <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = relhum_mean - 1.96 * relhum_sd,
                  ymax = relhum_mean + 1.96 * relhum_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = relhum_mean, color = "Measured microclimate", linetype = "Measured microclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = relhum_macro,
                color = "Measured macroclimate",
                linetype = "Measured macroclimate"),
            size = 0.6) +
  geom_line(aes(x = obs_time, y = relhum, color = "Simulated microclimate", linetype = "Simulated microclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = relhum_macro_sim,
                color = "CMIP6 macroclimate",
                linetype = "CMIP6 macroclimate"),
            size = 0.6) +
  scale_color_manual(name = "",
                     values = c("Measured microclimate" = "#3F826D",
                                "Measured macroclimate" = "#3F826D",
                                "Simulated microclimate" = "#FAC05E",
                                "CMIP6 macroclimate" = "#FAC05E")) +
  scale_linetype_manual(name = "",
                        values = c("Measured microclimate" = "solid",
                                   "Simulated microclimate" = "solid",
                                   "Measured macroclimate" = "dashed",
                                   "CMIP6 macroclimate" = "dashed")) +

  labs(x = "", y = "Relative Humidity (%)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 18),
    legend.title = element_text(face = "bold")
  )

# Print plots to a pdf file
pdf("../../figs/mc_output/mc_emp_vs_sim_mc_regua_era5_v1.pdf", height = 5, width = 10)
print((plot_airt + plot_relhum) +
  plot_layout(guides = "collect") &
  guides(color = guide_legend(ncol = 2)) &
  theme(legend.position = "bottom"))
dev.off()
