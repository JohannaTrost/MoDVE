library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)  # for smoothing later
library(viridis)  # good for perceptual color maps
library(scales)
library(lubridate)

indir = "/Users/johanna/Uni/masterarbeit/code/output/microclimf/all_days_sim_regua_2024"
mean_daily <- readRDS(paste(indir, "mean_daily_tz.rds", sep = "/"))
mean_daytime <- readRDS(paste(indir, "mean_daytime_tz.rds", sep = "/"))
mean_nighttime <- readRDS(paste(indir, "mean_nighttime_tz.rds", sep = "/"))

# Get ERA5 data
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))
# Calculate daily average temperature
daily_avg_temp <- climdata_reg %>%
  mutate(date = as.Date(obs_time)) %>%
  group_by(date) %>%
  summarize(daily_temp = mean(temp, na.rm = TRUE)) %>%
  arrange(date)
# Smooth
daily_avg_temp <- daily_avg_temp %>%
  mutate(daily_temp_smooth = zoo::rollmean(daily_temp, k = 10, fill = NA, align = "center"))

# Prepare data
daily_dates <- as.Date("2024-01-01") + 0:365
n_heights <- dim(mean_daily)[4]
heights <- seq(0.5, n_heights)
height_labels <- paste0(heights, " m")
height_colors <- viridis(n_heights)

# Extract data into long format dataframe
df <- data.frame(
  date = rep(daily_dates, times = n_heights),
  height = rep(heights, each = length(daily_dates)),
  Tz = as.vector(apply(mean_daily, c(3, 4), mean)),
  height_class = factor(rep(height_labels, each = length(daily_dates)))
)

# Apply 5-day moving average per height
df_smoothed <- df %>%
  group_by(height, height_class) %>%
  arrange(date) %>%
  mutate(Tz_smooth = zoo::rollmean(Tz, k = 10, fill = NA, align = "center"))

# Plot 2: Smoothed temperature curves
ggplot(df_smoothed, aes(x = date, y = Tz_smooth, color = height_class)) +
  geom_line() +
  scale_color_manual(
    values = setNames(height_colors, height_labels),
    name = "Height (m)"
  ) +
  labs(x = "Date", y = "Temperature") +
  theme_minimal() +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 5))) +
  # Add daily average temperature as red line:
  geom_line(data = daily_avg_temp, aes(x = date, y = daily_temp_smooth), color = "red", size = 0.5)

# Other plot 
ggplot(df_smoothed, aes(x = date, y = height, fill = Tz_smooth)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = median(df_smoothed$Tz_smooth, na.rm = TRUE),
    name = "Temperature"
  ) +
  labs(x = "Date", y = "Height (m)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


