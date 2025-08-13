
library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(lubridate)
library(furrr)
library(future)
library(microclimf)
library(dplyr)
library(ggplot2)
library(purrr)
library(readxl)
library(tidyverse)


extract_params <- function(raster_list, lon, lat, crs = "EPSG:4326") {
  # Create a SpatVector point in WGS84 (decimal degrees)
  point <- vect(cbind(lon, lat), crs = crs)

  # Transform the point to the CRS of the rasters
  target_crs <- crs(unwrap(raster_list[[1]]))
  point_proj <- project(point, target_crs)

  # Extract the value at the projected point for each raster
  result <- lapply(raster_list, function(r) {
    r <- unwrap(r)
    extract(r, point_proj)[[2]]
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

# Load data
region <- "pirineus"
outdir <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_output/", region, "_2024_test_v6")
if (region == "regua") {
  start_date <- "2024-10-27 00:00:00" # regua
  end_date <- "2024-10-31 23:59:59" # regua
} else if (region == "pirineus") {
  start_date <- "2024-09-20 00:00:00" # pirineus
  end_date <- "2024-09-23 23:59:59" # Pirineus
}

start <- 21
nx <- 4

sim_path <- paste0(outdir, "/sim_with_era5_",
                                        start + 1 , "_", start + nx, "_",
                                        substr(start_date, 1, 10),
                                        "_",
                                        substr(end_date, 1, 10),
                                        ".rds")
sim <- readRDS(sim_path)

# Dates
start_dt <- ymd_hms(start_date)
end_dt <- ymd_hms(end_date)

# Create hourly sequence as a vector
hourly_sequence <- seq(from = start_dt,
                       to = end_dt,
                       by = "hour")

# put into dataframe
mc_sim <- data.frame(
  obs_time = hourly_sequence,
  tair = apply(sim$tair[,, 2,], 3, mean),
  relhum = apply(sim$relhum[,, 2,], 3, mean)
)

# --- Load empirical data

if (region == "pirineus") {
  # --- PIRINEUS
  emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 770m elevation low Pirineus understorey transect"

  # Define logger directories
  emp_dirs <- c("10m", "30m", "50m", "70m", "90m")
  macro_dir <- "Reference or at 30m in rockfall sort of gap on pole"  # This is the macroclimate reference
} else if (region == "regua") {
  # --- REGUA
  emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 400m elevation REGUA understory Trilha Verde"
  emp_dirs <- c("3600m, 446mASL", "3800m, 387mASL", "3650m, 438mASL", "3700m, 433mASL")
  macro_dir <- "3850m, 387mASL waterfall reference"  # This is the macroclimate reference
}

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
      obs_hour_utc = ceiling_date(obs_time_utc, unit = "hour")
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
    mae_tair_model = mean(abs(joined$tair_emp - joined$tair), na.rm = TRUE),
    cor_tair_model = cor(joined$tair_emp, joined$tair, use = "complete.obs"),
    mae_relhum_model = mean(abs(joined$relhum_emp - joined$relhum), na.rm = TRUE),
    cor_relhum_model = cor(joined$relhum_emp, joined$relhum, use = "complete.obs"),

    # Macroclimate vs Empirical (skip if same as macro reference)
    mae_tair_macro = if (dir_name != macro_dir) mean(abs(joined$tair_emp - joined$tair_macro), na.rm = TRUE) else NA_real_,
    cor_tair_macro = if (dir_name != macro_dir) cor(joined$tair_emp, joined$tair_macro, use = "complete.obs") else NA_real_,
    mae_relhum_macro = if (dir_name != macro_dir) mean(abs(joined$relhum_emp - joined$relhum_macro), na.rm = TRUE) else NA_real_,
    cor_relhum_macro = if (dir_name != macro_dir) cor(joined$relhum_emp, joined$relhum_macro, use = "complete.obs") else NA_real_
  )
}

# Apply to all loggers
results <- map_dfr(emp_dirs, process_comparison)

# View sorted by model RH MAE
results_sorted <- results %>% arrange(mae_relhum_model)
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
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) #%>%
  #rename(!!repl_tair := "tair", !!repl_relhum := "relhum")

names(emp_data_combined)

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
  inner_join(., mc_sim, by = "obs_time")

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = tair_mean - 1.96 * tair_sd,
                  ymax = tair_mean + 1.96 * tair_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = tair_mean, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = tair_macro, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = tair, color = "Simulated"),
            size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01"),       # Orange
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated")) +
  labs(title = "Air Temperature: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Air Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Second plot: Relative Humidity (empirical vs. simulated)
plot_relhum <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = relhum_mean - 1.96 * relhum_sd,
                  ymax = relhum_mean + 1.96 * relhum_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = relhum_mean, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = relhum_macro, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = relhum, color = "Simulated"),
            size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01"),       # Orange
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated")) +
  labs(title = "Relative Humidity: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Relative Humidity (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Print plots to a pdf file
pdf(paste0("../../figs/mc_output/airt_emp_vs_sim_mc_", region, "_v6.pdf"))
print(plot_airt)
dev.off()

pdf(paste0("../../figs/mc_output/relhum_emp_vs_sim_mc_", region, "_v6.pdf"))
print(plot_relhum)
dev.off()


# pdf(paste0("../../figs/mc_output/relhum_emp_vs_sim_mc_", region, "_ccf_v6.pdf"))
# ccf(emp_sim_data$tair, emp_sim_data$tair_emp)
# dev.off()

# ----- Compare on daily scale

# Aggregate hourly data to daily by averaging
emp_daily <- emp_summary %>%
  mutate(date = as_date(obs_time)) %>%  # Extract date from datetime
  group_by(date) %>%
  summarise(
    tair_mean_daily = mean(tair_mean, na.rm = TRUE),
    tair_sd_daily = mean(tair_sd, na.rm = TRUE),
    relhum_mean_daily = mean(relhum_mean, na.rm = TRUE),
    relhum_sd_daily = mean(relhum_sd, na.rm = TRUE),
    tair_macro_daily = mean(tair_macro, na.rm = TRUE),
    relhum_macro_daily = mean(relhum_macro, na.rm = TRUE),
    tair_daily = mean(tair, na.rm = TRUE),
    relhum_daily = mean(relhum, na.rm = TRUE),
    n_hours = n()  # Number of hourly observations per day
  ) %>%
  ungroup()

# View the result
print(emp_daily)

# Add CMIP and ERA5 data to daily summary
indir <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_input/", region)
cmip_path <- paste(in_dir, "cmip6_climdata_2024_v2.csv", sep = "/")
cmip <- read_csv(cmip_path)

era5_path <- paste(in_dir, "era5_climdata_2024.csv", sep = "/")
era5 <- read_csv(era5_path)

# Aggregate CMIP and ERA5 to daily
cmip_daily <- cmip %>%
  mutate(date = as_date(obs_time)) %>%  # Extract date from datetime
  group_by(date) %>%
  summarise(
    tair_cmip = mean(temp, na.rm = TRUE),
    relhum_cmip = mean(relhum, na.rm = TRUE),
  ) %>%
  ungroup()

era5_daily <- era5 %>%
  mutate(date = as_date(obs_time)) %>%  # Extract date from datetime
  group_by(date) %>%
  summarise(
    tair_era5 = mean(temp, na.rm = TRUE),
    relhum_era5 = mean(relhum, na.rm = TRUE),
  ) %>%
  ungroup()

# Join CMIP and ERA5 data with empirical daily data
all_daily <- emp_daily %>%
  left_join(cmip_daily, by = "date") %>%
  left_join(era5_daily, by = "date")

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(all_daily) +
  geom_ribbon(aes(x = date,
                  ymin = tair_mean_daily - 1.96 * tair_sd_daily,
                  ymax = tair_mean_daily + 1.96 * tair_sd_daily),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = date, y = tair_mean_daily, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = date, y = tair_macro_daily, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = date, y = tair_cmip, color = "CMIP6 Ensemble"),
            size = 1) +
  geom_line(aes(x = date, y = tair_era5, color = "ERA5"), size = 1) +
  geom_line(aes(x = date, y = tair_daily, color = "Simulated"),
            size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01",       # Orange
                                "CMIP6 Ensemble" = "#1F77B4",   # Blue
                                "ERA5" = "red"),            # red
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated",
                                "CMIP6 Ensemble", "ERA5")) +
  labs(title = "Air Temperature: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Air Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Second plot: Relative Humidity (empirical vs. simulated)
plot_relhum <- ggplot(all_daily) +
  geom_ribbon(aes(x = date,
                  ymin = relhum_mean_daily - 1.96 * relhum_sd_daily,
                  ymax = relhum_mean_daily + 1.96 * relhum_sd_daily),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = date, y = relhum_mean_daily, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = date, y = relhum_macro_daily, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = date, y = relhum_daily, color = "Simulated"),
            size = 1) +
  geom_line(aes(x = date, y = relhum_cmip, color = "CMIP6 Ensemble"),
            size = 1) +
  geom_line(aes(x = date, y = relhum_era5, color = "ERA5"), size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01",       # Orange
                                "CMIP6 Ensemble" = "#1F77B4",   # Blue
                                "ERA5" = "red"),            # red
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated",
                                "CMIP6 Ensemble", "ERA5")) +
  labs(title = "Relative Humidity: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Relative Humidity (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Print plots to a pdf file
pdf(paste0("../../figs/mc_output/airt_daily_emp_vs_sim_mc_", region, "_v6.pdf"))
print(plot_airt)
dev.off()

pdf(paste0("../../figs/mc_output/relhum_daily_emp_vs_sim_mc_", region, "_v6.pdf"))
print(plot_relhum)
dev.off()

# Final stats

all_daily %>%
  summarise(
    tair_MacroObs = mean(abs(tair_macro_daily - tair_daily), na.rm = TRUE),
    tair_CMIP  = mean(abs(tair_cmip - tair_daily), na.rm = TRUE),
    tair_ERA5  = mean(abs(tair_era5 - tair_daily), na.rm = TRUE),
    tair_MicroObs = mean(abs(tair_mean_daily - tair_daily), na.rm = TRUE),

    relhum_MacroObs = mean(abs(relhum_macro_daily - relhum_daily), na.rm = TRUE),
    relhum_CMIP  = mean(abs(relhum_cmip - relhum_daily), na.rm = TRUE),
    relhum_ERA5  = mean(abs(relhum_era5 - relhum_daily), na.rm = TRUE),
    relhum_MicroObs = mean(abs(relhum_mean_daily - relhum_daily), na.rm = TRUE)
  ) %>%
pivot_longer(
    cols = everything(),
    names_to = "Comparison",
    values_to = "MAE"
  )