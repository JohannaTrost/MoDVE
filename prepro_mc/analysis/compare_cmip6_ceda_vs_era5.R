library(readr)      # for reading NetCDF
library(dplyr)      # tidyverse for data wrangling
library(lubridate)  # for datetime handling
library(purrr)      # for map functions
library(tidyr)
library(terra)
library(ggplot2)

cmip_dir <- file.path("../../data/mc_input/climate/cmip6_ceda")
# Load Cmipt climate data
cmip6_temp_ts <- read_csv(file.path(cmip_dir, "cmip6_ceda_tas_1981_2024.csv"))
cmip6_rh_ts   <- read_csv(file.path(cmip_dir, "cmip6_ceda_hurs_1981_2024.csv"))

# Load ERA5 climate data
in_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/mc_input/climate/era5_processed")

era5_temp_ts <- tibble()
era5_rh_ts <- tibble()

for (year in c(1981:2024)) {
  era5_raw <- readRDS(file.path(in_dir, paste0("era5_climdata_", year, ".RDS")))
  era5_temp <- terra::unwrap(era5_raw$temp)
  era5_relhum <- terra::unwrap(era5_raw$relhum)

  layer_names <- names(terra::unwrap(era5_raw$precip))
  timestamps <- sub(".*=", "", layer_names)
  timestamps <- as.POSIXct(as.numeric(timestamps), origin = "1970-01-01", tz = "UTC")

  era5_avg_tas <- tibble(
    time = timestamps,
    tas = terra::global(era5_temp, "mean", na.rm = TRUE)[,1],
    model = "ERA5",
  )
  era5_avg_rh <- tibble(
    time = timestamps,
    hurs = terra::global(era5_relhum, "mean", na.rm = TRUE)[,1],
    model = "ERA5",
  )
  era5_temp_ts <- bind_rows(era5_temp_ts, era5_avg_tas)
  era5_rh_ts <- bind_rows(era5_rh_ts, era5_avg_rh)
}

# Aggregate ERA5 from hourly to daily
era5_temp_daily <- era5_temp_ts %>%
  mutate(time = as.Date(time)) %>%
  group_by(time, model) %>%
  summarize(tas = mean(tas, na.rm = TRUE), .groups = "drop")
era5_rh_daily <- era5_rh_ts %>%
    mutate(time = as.Date(time)) %>%
    group_by(time, model) %>%
    summarize(hurs = mean(hurs, na.rm = TRUE), .groups = "drop")

# Concatenaete CMIP6 and ERA5 data
temp_ts <- bind_rows(
  cmip6_temp_ts,
  era5_temp_daily
)
rh_ts <- bind_rows(
  cmip6_rh_ts,
  era5_rh_daily
)

# Aggregate combined data to monthly
temp_monthly <- temp_ts %>%
  mutate(year = year(time), month = month(time)) %>%
  group_by(year, month, model) %>%
  summarize(tas = mean(tas, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste(year, month, "15", sep = "-"))) %>%
  select(date, year, month, model, tas)

# For debugging print the count of each month for all years for each model
print(temp_monthly %>%
  group_by(year) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(year) %>%
  print(n = Inf))

# Aggregate combined data to monthly fro rh
rh_monthly <- rh_ts %>%
  mutate(year = year(time), month = month(time)) %>%
  group_by(year, month, model) %>%
  summarize(hurs = mean(hurs, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste(year, month, "15", sep = "-"))) %>%
  select(date, year, month, model, hurs)

# ------------ PLOT RH MONTHLY ------------

# Filter data for ERA5 and CMIP6 Average models only
plot_data <- rh_monthly %>%
  filter(model %in% c("ERA5", "CMIP6 Average"))

# Create color palettes
# ERA5: Light green (1981) to dark green (2024)
era5_colors <- colorRampPalette(c("#90EE90", "#006400"))(length(unique(plot_data$year)))

# CMIP6: Dark blue (1981) to light blue (2024)
cmip6_colors <- colorRampPalette(c("#00008B", "#87CEEB"))(length(unique(plot_data$year)))

# Create the plot
p <- ggplot(plot_data, aes(x = month, y = hurs, color = interaction(model, year))) +
  geom_line(aes(group = interaction(model, year)), alpha = 0.7, size = 0.8) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  ) +
  labs(
    title = "Monthly Relative Humidity by Model and Year (1981-2024)",
    x = "Month",
    y = "Relative Humidity (%)",
    caption = "Light to dark green: ERA5 (1981-2024)\nDark to light blue: CMIP6 Average (1981-2024)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",  # Remove legend since we have many lines
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid.minor = element_blank()
  )

# Manually assign colors based on model and year
years <- unique(plot_data$year)
color_mapping <- c()

for (i in seq_along(years)) {
  year <- years[i]
  color_mapping[paste0("ERA5.", year)] <- era5_colors[i]
  color_mapping[paste0("CMIP6 Average.", year)] <- cmip6_colors[i]
}

p <- p + scale_color_manual(values = color_mapping)

# Print both versions
pdf("../../figs/mc_input/compare_rh_cmip6_era5_monthly_1981_2024_v1.pdf")
print(p)
dev.off()

# ------------ PLOT TEMP MONTHLY ------------

# Filter data for ERA5 and CMIP6 Average models only
plot_data <- temp_monthly %>%
  filter(model %in% c("ERA5", "CMIP6 Average"))

# Create the plot
p <- ggplot(plot_data, aes(x = month, y = tas, color = interaction(model, year))) +
  geom_line(aes(group = interaction(model, year)), alpha = 0.7, size = 0.8) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  ) +
  labs(
    title = "Monthly Temperature by Model and Year (1981-2024)",
    x = "Month",
    y = "Temperature (°C)",
    caption = "Light to dark green: ERA5 (1981-2024)\nDark to light blue: CMIP6 Average (1981-2024)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",  # Remove legend since we have many lines
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid.minor = element_blank()
  )

# Manually assign colors based on model and year
years <- unique(plot_data$year)
color_mapping <- c()

for (i in seq_along(years)) {
  year <- years[i]
  color_mapping[paste0("ERA5.", year)] <- era5_colors[i]
  color_mapping[paste0("CMIP6 Average.", year)] <- cmip6_colors[i]
}

p <- p + scale_color_manual(values = color_mapping)

# Print both versions
pdf("../../figs/mc_input/compare_temp_cmip6_era5_monthly_1981_2024_v1.pdf")
print(p)
dev.off()

# ------------------ Yealy aggregates and plots ------------------

temp_annual <- temp_monthly %>%
  group_by(year, model) %>%
  summarise(tas = mean(tas, na.rm = TRUE), .groups = "drop")

rh_annual <- rh_monthly %>%
  group_by(year, model) %>%
  summarise(hurs = mean(hurs, na.rm = TRUE), .groups = "drop")

# Plot annual temperature

# Calculate standard deviation for CMIP6 models (excluding ERA5 and CMIP6 Average)
cmip6_models_sd <- temp_annual %>%
  filter(!model %in% c("ERA5", "CMIP6 Average")) %>%
  group_by(year) %>%
  summarise(sd_tas = sd(tas, na.rm = TRUE), .groups = 'drop')

# Get ERA5 and CMIP6 Average data
main_data <- temp_annual %>%
  filter(model %in% c("ERA5", "CMIP6 Average"))

# Merge with standard deviation data
plot_data <- main_data %>%
  left_join(cmip6_models_sd, by = "year") %>%
  mutate(
    upper_bound = ifelse(model == "CMIP6 Average", tas + sd_tas, NA),
    lower_bound = ifelse(model == "CMIP6 Average", tas - sd_tas, NA)
  )

# Create the plot
p <- ggplot(plot_data, aes(x = year, y = tas)) +
  # Add uncertainty ribbon for CMIP6 Average
  geom_ribbon(data = filter(plot_data, model == "CMIP6 Average"),
              aes(ymin = lower_bound, ymax = upper_bound),
              fill = "lightblue", alpha = 0.3) +

  # Add lines for both models
  geom_line(aes(color = model, linetype = model), size = 1) +

  # Add points
  geom_point(aes(color = model, shape = model), size = 2) +

  # Customize colors and shapes
  scale_color_manual(values = c("ERA5" = "forestgreen", "CMIP6 Average" = "steelblue")) +
  scale_linetype_manual(values = c("ERA5" = "solid", "CMIP6 Average" = "solid")) +
  scale_shape_manual(values = c("ERA5" = 16, "CMIP6 Average" = 17)) +

  # Labels and theme
  labs(
    title = "Annual Temperature: ERA5 vs CMIP6 Average with Model Uncertainty",
    x = "Year",
    y = "Temperature (°C)",
    color = "Model",
    linetype = "Model",
    shape = "Model",
    caption = "Shaded area represents ±1 standard deviation of CMIP6 models"
  ) +

  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Display the plot
pdf("../../figs/mc_input/compare_temp_cmip6_era5_annual_1981_2024_v1.pdf")
print(p)
dev.off()

# Plot annual relative humidity

# Calculate standard deviation for CMIP6 models (excluding ERA5 and CMIP6 Average)
cmip6_models_sd <- rh_annual %>%
  filter(!model %in% c("ERA5", "CMIP6 Average")) %>%
  group_by(year) %>%
  summarise(sd_hurs = sd(hurs, na.rm = TRUE), .groups = 'drop')

# Get ERA5 and CMIP6 Average data
main_data <- rh_annual %>%
  filter(model %in% c("ERA5", "CMIP6 Average"))

# Merge with standard deviation data
plot_data <- main_data %>%
  left_join(cmip6_models_sd, by = "year") %>%
  mutate(
    upper_bound = ifelse(model == "CMIP6 Average", hurs + sd_hurs, NA),
    lower_bound = ifelse(model == "CMIP6 Average", hurs - sd_hurs, NA)
  )

# Create the plot
p <- ggplot(plot_data, aes(x = year, y = hurs)) +
  # Add uncertainty ribbon for CMIP6 Average
  geom_ribbon(data = filter(plot_data, model == "CMIP6 Average"),
              aes(ymin = lower_bound, ymax = upper_bound),
              fill = "lightblue", alpha = 0.3) +

  # Add lines for both models
  geom_line(aes(color = model, linetype = model), size = 1) +

  # Add points
  geom_point(aes(color = model, shape = model), size = 2) +

  # Customize colors and shapes
  scale_color_manual(values = c("ERA5" = "forestgreen", "CMIP6 Average" = "steelblue")) +
  scale_linetype_manual(values = c("ERA5" = "solid", "CMIP6 Average" = "solid")) +
  scale_shape_manual(values = c("ERA5" = 16, "CMIP6 Average" = 17)) +

  # Labels and theme
  labs(
    title = "Annual Relative Humidity: ERA5 vs CMIP6 Average with Model Uncertainty",
    x = "Year",
    y = "Relative Humidity (%)",
    color = "Model",
    linetype = "Model",
    shape = "Model",
    caption = "Shaded area represents ±1 standard deviation of CMIP6 models"
  ) +

  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Display the plot
pdf("../../figs/mc_input/compare_rh_cmip6_era5_annual_1981_2024_v1.pdf")
print(p)
dev.off()
