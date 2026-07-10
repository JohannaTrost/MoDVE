# ------
# Validate microclimate using measurements at Johansson Zones in REGUA, Atlantic forest (2025)

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
library(DEoptim)
library(parallel)
library(doFuture)
library(progressr)
library(foreach)

# CONFIGURE paths
emp_dir <- file.path("../modve_data_zenodo/empirical/300_500 REGUA")
figs_out <- file.path("../modve_figs/mc_output")
in_dir <- file.path("..", "modve_data_zenodo", "mc_input")
in_dir_regua <- file.path("..", "modve_data_zenodo", "mc_input", "regua")
veg_path <- file.path(in_dir_regua, "rep0", paste0("vegp_mof3d_ptm_124.RDS"))
soil_path <- paste(in_dir_regua, "soil.RDS", sep = "/")
clim_path <- file.path(in_dir, "climate", "era5_processed", "era5_climdata_2025.csv") # ERA5
mh_dir <- file.path("../modve_data_zenodo/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/microhabitat_mc") # Regua

format_time <- function(secs) {
  secs <- max(0L, round(secs))
  if (secs < 60)   return(sprintf("%ds", secs))
  if (secs < 3600) return(sprintf("%dm %02ds", secs %/% 60L, secs %% 60L))
  sprintf("%dh %02dm", secs %/% 3600L, (secs %% 3600L) %/% 60L)
}


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

##################################################################################
#                               Empirical data                                   #
##################################################################################

# ---- Heights of loggers

# Read the file (assuming it's named "data.txt" in your working directory)
heights_measured <- read_lines(file.path(emp_dir, "300_500_REGUA_survey_heights.txt"))

# Extract the relevant lines (JZ1, JZ2, etc.)
logger_data <- heights_measured[str_detect(heights_measured, "JZ")]

# Split into logger and height, then clean
heights_measured <- tibble(
  line = logger_data
) %>%
  separate(line, into = c("logger", "height"), sep = " ", remove = FALSE) %>%
  mutate(
    height = str_remove(height, "m"),  # Remove 'm'
    height = as.numeric(height)         # Convert to numeric
  ) %>%
  select(logger, height)               # Keep only the two columns

# View the result
print(heights_measured)

# ---- Read in logger data

# Define logger directories
emp_dirs <- c("JZ1", "JZ2", "JZ3", "JZ4", "JZ5")
macro_dir <- "Control"  # This is the macroclimate reference

# Function to load and process a logger file
load_logger_data <- function(dir_name) {
  file_path <- file.path(emp_dir, dir_name, "mc_data.xlsx")
  file_path_temp <- file.path(emp_dir, dir_name, "temp.xlsx")

  temp_df <- read_excel(file_path_temp) %>%  select("Date-Time (Brazil Standard Time)", "Temperature , °C")
  light_df <- read_excel(file_path_temp) %>%  select("Date-Time (Brazil Standard Time)", "Light , lux")
  data <- left_join(temp_df, light_df, by = "Date-Time (Brazil Standard Time)")

  if (file.exists(file_path)) {
    hum_df <- read_excel(file_path) %>% select("Date-Time (Brazil Standard Time)", "RH , %")
    data <- left_join(data, hum_df, by = "Date-Time (Brazil Standard Time)")
  } else {
    data["RH , %"] <- NA
  }

  data %>%
    select("Date-Time (Brazil Standard Time)", "Light , lux", "Temperature , °C", "RH , %") %>%
    rename(obs_time = "Date-Time (Brazil Standard Time)",
           tair = "Temperature , °C",
           relhum = "RH , %",
           light = "Light , lux"
    ) %>%
    mutate(
      obs_time = force_tz(as.POSIXct(obs_time), tzone = "America/Sao_Paulo"),
      obs_time_utc = with_tz(obs_time, tzone = "UTC"),
      obs_hour_utc = floor_date(obs_time_utc, unit = "hour")
    ) %>%
    group_by(obs_hour_utc) %>%
    summarise(
      tair = mean(tair, na.rm = TRUE),
      relhum = mean(relhum, na.rm = TRUE),
      light = mean(light, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(obs_time = obs_hour_utc)
}

# Load macroclimate data
macro_data <- load_logger_data(macro_dir) %>%
  rename(tair_macro = tair, relhum_macro = relhum, light_macro = light)

# -- Plot emp. data

# Load and process all logger data
emp_data_list <- map(emp_dirs, load_logger_data)

# Combine all into a single tibble by full_join by time
repl_tair <- paste0("tair_", emp_dirs[length(emp_dirs)])
repl_relhum <- paste0("relhum_", emp_dirs[length(emp_dirs)])
repl_light <- paste0("light_", emp_dirs[length(emp_dirs)])

emp_data_combined <- emp_data_list %>%
  reduce(full_join, by = "obs_time") %>%
  rename_at(vars(starts_with("tair")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename_at(vars(starts_with("relhum")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename_at(vars(starts_with("light")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename(!!repl_tair := "tair", !!repl_relhum := "relhum", !!repl_light := "light")

names(emp_data_combined)

# Reshape data to long format for both temperature and humidity
emp_data_long <- emp_data_combined %>%
  pivot_longer(
    cols = c(tair_JZ1, tair_JZ2, tair_JZ3, tair_JZ4, tair_JZ5,
             relhum_JZ1, relhum_JZ2, relhum_JZ3, relhum_JZ4, relhum_JZ5,
             light_JZ1, light_JZ2, light_JZ3, light_JZ4, light_JZ5
    ),
    names_to = c("variable", "height"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    variable = case_when(
      variable == "tair" ~ "Temperature (°C)",
      variable == "relhum" ~ "Relative Humidity (%)",
      variable == "light" ~ "Light (lux)",
      TRUE ~ variable
    ),
    height = case_when(
      height == "JZ1" ~ "0.4m (JZ1)",
      height == "JZ2" ~ "5.5m (JZ2)",
      height == "JZ3" ~ "9.5m (JZ3)",
      height == "JZ4" ~ "12.5m (JZ4)",
      height == "JZ5" ~ "15.5m (JZ5)",
      TRUE ~ height
    )
  )

# Define colors for each height
height_colors <- c(
  "0.4m (JZ1)" = "dodgerblue",
  "5.5m (JZ2)" = "firebrick",
  "9.5m (JZ3)" = "darkorange",
  "12.5m (JZ4)" = "forestgreen",
  "15.5m (JZ5)" = "purple"
)

# Plot Microclimate time series at all JZ
combined_plot <- ggplot(emp_data_long %>% filter(variable != "Light (lux)"),
                        aes(x = obs_time, y = value, color = height, linetype = height)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = height_colors) +
  scale_linetype_manual(values = rep(1, 5)) +  # Solid lines for all
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  labs(
    x = "",
    y = "",
    color = "Height",
    linetype = "Height"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 18),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )

# Print the combined plot
pdf(file.path(figs_out, "2025_empirical_jz_mc_2025_regua.pdf"), height = 5, width = 10)
print(combined_plot)
dev.off()

# For later comparisons
emp_data <- emp_data_long %>%
  # Pivot to wider format
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  # Rename columns for clarity
  rename("tair_emp" = "Temperature (°C)", "relhum_emp" = "Relative Humidity (%)", "light_emp" = "Light (lux)") %>%
  mutate(logger = sub(".*\\(([A-Za-z0-9]+)\\)", "\\1", height)) %>%
  select(-height)

##################################################################################
#                               Simulated data                                   #
##################################################################################

# Set Parameters
optimize <- FALSE
x <- 25
y <- 25
ts <- 123
forst <- 0

# Load data
vegp_reg <- readRDS(veg_path)  # PAI and canopy height will be replaced by MoF3D output below
soilc_reg <- readRDS(soil_path)
climdata_reg <- read_csv(clim_path) # REGUA

# Get coordiantes
coords_veg <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
coords_soil <- indices2coords(x, y, terra::unwrap(soilc_reg$aspect))[c("x", "y")]
lon <- coords_soil[[1]]
lat <- coords_soil[[2]]

# Get params for the point
vegparams <- extract_params(vegp_reg, coords_veg[[1]], coords_veg[[2]])
grndparams <- extract_params(soilc_reg, lon, lat)

# Update veg. height to actual height
vegparams$h <- 22

# Get PAI
microhab_file <- file.path(mh_dir, paste0("rep", forst), paste0("microhabitatMatrix", ts, ".rds")) # Regua
pai <- readRDS(microhab_file)[,,,1]
paii <- apply(pai[,,1:max(vegparams$h, 0.5)], c(3), mean, na.rm = TRUE)
vegparams$pai <- sum(paii)

# --- Run the simulation with the paii

# Run simulations for all heights/loggers with optimized paii
mc_sim <- data.frame()
for (i in seq_along(heights_measured$height)) {
  h <- heights_measured$height[i]
  jz <- heights_measured$logger[i]

  mout <- micropoint::runpointmodel(
    climdata_reg, reqhgt = h, vegparams,
    paii, grndparams, lat = lat, long = lon
  )

  mc_sim_h <- data.frame(
    obs_time = mout$obs_time,
    tair = mout$tair,
    tleaf = mout$tleaf,
    relhum = mout$relhum,
    logger = jz,
    height = h
  )

  # Filter to empirical dates
  mc_sim_h_filtered <- mc_sim_h[mc_sim_h$obs_time %in% macro_data$obs_time, ]
  mc_sim <- rbind(mc_sim, mc_sim_h_filtered)
}

# Join with empirical and macro data
model_macro <- climdata_reg %>%
  select(obs_time, temp, relhum) %>%
  rename(tair_macro_sim = temp, relhum_macro_sim = relhum)

joined <- emp_data %>%
  left_join(mc_sim, by = c("obs_time", "logger")) %>%
  left_join(macro_data, by = "obs_time") %>%
  left_join(model_macro, by = "obs_time") %>%
  rename("JZ" = "logger")

safe_cor_test <- function(x, y) {
  if (all(is.na(x)) | all(is.na(y))) {
    return(list(cor = NA, p.value = NA))
  } else {
    test <- cor.test(x, y, use = "complete.obs")
    return(list(cor = test$estimate, p.value = test$p.value))
  }
}

# Compute RMSE for each JZ/height
joined %>%
  group_by(JZ, height) %>%
  summarize(
    rmse_tair = sqrt(mean((tair_emp - tair)^2, na.rm = TRUE)),
    rmse_relhum = sqrt(mean((relhum_emp - relhum)^2, na.rm = TRUE)),
    mae_tair = mean(tair_emp - tair, na.rm = TRUE),
    mae_relhum = mean(relhum_emp - relhum, na.rm = TRUE),
    cor_tair = safe_cor_test(tair_emp, tair)$cor,
    cor_relhum = safe_cor_test(relhum_emp, relhum)$cor,
    .groups = "drop"
  )

# ------------------ Plot understorey/upper canopy ts only

h <- 0.4 # also do for 15.5

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(joined %>% filter(height == h)) +
  geom_line(aes(x = obs_time, y = tair, color = "Simulated microclimate", linetype = "Simulated microclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = tair_emp, color = "Measured microclimate", linetype = "Measured microclimate"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = tair_macro,
                color = "Measured macroclimate",
                linetype = "Measured macroclimate"),
            size = 0.6) +
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
plot_relhum <- ggplot(joined %>% filter(height == h)) +
  geom_line(aes(x = obs_time, y = relhum_emp, color = "Measured microclimate", linetype = "Measured microclimate"),
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
pdf(file.path(figs_out, paste0("mc_emp_vs_sim_regua_era5_2025_", h, "m.pdf")), height = 5, width = 10)
print((plot_airt + plot_relhum) +
 plot_layout(guides = "collect") &
 guides(color = guide_legend(ncol = 2)) &
 theme(legend.position = "bottom"))
dev.off()

# ------------------ Plot profile

# Average over days for 3 specific points in time (8am, 12am, 6pm)
gradient_agg <- joined %>%
  filter(!date(obs_time) %in% as.Date(c("2025-09-25", "2025-09-26"))) %>%
  mutate(time = hour(obs_time)) %>%
  group_by(height, time) %>%
  summarize(
    relhum_emp_mean = mean(relhum_emp, na.rm = TRUE),
    relhum_mean = mean(relhum, na.rm = TRUE),
    relhum_macro_mean = mean(relhum_macro, na.rm = TRUE),
    relhum_macro_sim_mean = mean(relhum_macro_sim, na.rm = TRUE),
    tair_emp_mean = mean(tair_emp, na.rm = TRUE),
    tair_mean = mean(tair, na.rm = TRUE),
    tair_macro_mean = mean(tair_macro, na.rm = TRUE),
    tair_macro_sim_mean = mean(tair_macro_sim, na.rm = TRUE),
    relhum_emp_sd = sd(relhum_emp, na.rm = TRUE),
    relhum_sd = sd(relhum, na.rm = TRUE),
    relhum_macro_sd = sd(relhum_macro, na.rm = TRUE),
    relhum_macro_sim_sd = sd(relhum_macro_sim, na.rm = TRUE),
    tair_emp_sd = sd(tair_emp, na.rm = TRUE),
    tair_sd = sd(tair, na.rm = TRUE),
    tair_macro_sd = sd(tair_macro, na.rm = TRUE),
    tair_macro_sim_sd = sd(tair_macro_sim, na.rm = TRUE),
  )

# ---- Plot gradient for each hour across all days

# Define hour order starting at 8am
hour_order <- c(8:23, 0:7)

# Prepare data: filter to hours with data and apply ordered factor
gradient_facet <- gradient_agg %>%
  filter(time %in% hour_order) %>%
  mutate(time_label = factor(paste0(sprintf("%02d", time), ":00"),
                             levels = paste0(sprintf("%02d", hour_order), ":00")))

grad_temp_plt_all_hours <- ggplot(gradient_facet, aes(y = height)) +
  geom_path(aes(x = tair_mean, color = "Simulated microclimate"), size = 1) +
  geom_ribbon(aes(xmin = tair_mean - 1.96 * tair_sd,
                  xmax = tair_mean + 1.96 * tair_sd,
                  fill = "±Std.Err.(simulated microclimate)"), alpha = 0.2) +
  geom_path(aes(x = tair_emp_mean, color = "Measured microclimate"), size = 1) +
  geom_ribbon(aes(xmin = tair_emp_mean - 1.96 * tair_emp_sd,
                  xmax = tair_emp_mean + 1.96 * tair_emp_sd,
                  fill = "±Std.Err.(measured microclimate)"), alpha = 0.2) +

  # Measured macroclimate: dashed boundary lines instead of ribbon
  geom_path(aes(x = tair_macro_mean, color = "Measured macroclimate"), size = 1, linetype = "dashed") +
  geom_path(aes(x = tair_macro_mean - 1.96 * tair_macro_sd, color = "Measured macroclimate"),
            size = 0.5, linetype = "dotted") +
  geom_path(aes(x = tair_macro_mean + 1.96 * tair_macro_sd, color = "Measured macroclimate"),
            size = 0.5, linetype = "dotted") +

  # CMIP6 macroclimate: dashed boundary lines instead of ribbon
  geom_path(aes(x = tair_macro_sim_mean, color = "CMIP6 macroclimate"), size = 1, linetype = "dashed") +
  geom_path(aes(x = tair_macro_sim_mean - 1.96 * tair_macro_sim_sd, color = "CMIP6 macroclimate"),
            size = 0.5, linetype = "dotted") +
  geom_path(aes(x = tair_macro_sim_mean + 1.96 * tair_macro_sim_sd, color = "CMIP6 macroclimate"),
            size = 0.5, linetype = "dotted") +

  scale_color_manual(
    name = "Temperature Gradient",
    values = c("Measured microclimate"  = "#3F826D",
               "Measured macroclimate"  = "#3F826D",
               "Simulated microclimate" = "#FAC05E",
               "CMIP6 macroclimate"     = "#FAC05E")
  ) +
  scale_fill_manual(
    name = "Std. Err.",
    values = c("±Std.Err.(measured microclimate)"  = "#3F826D",
               "±Std.Err.(simulated microclimate)" = "#FAC05E")
  ) +
  facet_wrap(~ time_label, ncol = 6) +
  labs(
    x = "Temperature",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 12)
  )


pdf(file.path(figs_out, paste0("mc_emp_vs_sim_regua_era5_2025_gradient_hourly.pdf")), height = 16, width = 20)
print(grad_temp_plt_all_hours)
dev.off()

# ---- Plot gradient for specific day

for (day in c(26, 27, 28)) {
  # Filter for specific day
  gradient_single_day <- joined %>%
    filter(date(obs_time) == as.Date(paste0("2025-09-", day))) %>%
    mutate(time_label = factor(paste0(sprintf("%02d", hour(obs_time)), ":00"),
                               levels = paste0(sprintf("%02d", hour_order), ":00")))

  grad_temp_plt_single_day <- ggplot(gradient_single_day, aes(y = height)) +
    geom_path(aes(x = tair, color = "Simulated microclimate"), size = 1) +
    geom_path(aes(x = tair_emp, color = "Measured microclimate"), size = 1) +
    geom_path(aes(x = tair_macro, color = "Measured macroclimate"), size = 1, linetype = "dashed") +
    geom_path(aes(x = tair_macro_sim, color = "ERA5 macroclimate"), size = 1, linetype = "dashed") +
    scale_color_manual(
      name = "Temperature Gradient",
      values = c("Measured microclimate"  = "#3F826D",
                 "Measured macroclimate"  = "#3F826D",
                 "Simulated microclimate" = "#FAC05E",
                 "ERA5 macroclimate"     = "#FAC05E")
    ) +
    facet_wrap(~ time_label, ncol = 6) +
    labs(
      x = "Temperature (°C)",
      y = "Height (m)",
      title = paste0(day, " Sep 2025")
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 24),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 28),
      plot.title = element_text(size = 28),  # Double the title size
      legend.text = element_text(size = 24),  # Double the legend text size
      legend.title = element_blank()  # Double the legend title size
    )


  pdf(file.path(figs_out, paste0("mc_emp_vs_sim_regua_era5_2025_gradient_", day, ".pdf")), height = 16, width = 16)
  print(grad_temp_plt_single_day)
  dev.off()
}
