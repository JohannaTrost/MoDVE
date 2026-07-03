# Hourly Linear Detrending for Temperature Data
# This code removes long-term climate trends while preserving diurnal and seasonal cycles

library(dplyr)
library(lubridate)
library(readr)

# --- Scenario: 119 years until 2100 without CC - REGUA

region <- "regua"
data_dir <- file.path("../modve_data")
figs_dir <- file.path("..", "modve_figs", "mc_input")
out_dir <- file.path(data_dir, "mc_input", region, "scenarios")
climdata <- read_csv(file.path(out_dir, "climdata_era5_cmip6_1981-2100_ssp245.csv"))

if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
}

# Function to perform hourly linear detrending
detrend_hourly_mc <- function(data, mc_col = "temp", time_col = "obs_time",
                                       baseline_years = 1981:1985) {

  # Create additional time variables for grouping
  data_with_time <- data %>%
    mutate(
      year = year(!!sym(time_col)),
      month = month(!!sym(time_col)),
      day = day(!!sym(time_col)),
      hour = hour(!!sym(time_col)),
      # Create day of year (handling leap years consistently)
      day_of_year = yday(!!sym(time_col)),
      # Create unique identifier for each hour of each day of year
      day_hour = paste(day_of_year, sprintf("%02d", hour), sep = "_")
    )

  # Perform detrending for each unique day-hour combination
  detrended_data <- data_with_time %>%
    group_by(day_hour, month, day, hour) %>%
    mutate(
      detrended = {
        # Skip if insufficient data points (need at least 3 points for meaningful trend)
        if (n() < 3 || sum(!is.na(!!sym(mc_col))) < 3) {
          !!sym(mc_col)
        } else {
          # Remove missing values for trend fitting
          valid_indices <- !is.na(!!sym(mc_col))

          if (sum(valid_indices) < 3) {
            !!sym(mc_col)
          } else {
            # Fit linear trend: temperature ~ year
            values <- !!sym(mc_col)
            year_values <- year

            # Fit trend only on valid data
            valid_temp <- values[valid_indices]
            valid_years <- year_values[valid_indices]

            trend_model <- lm(valid_temp ~ valid_years)

            # Calculate trend values for all data points
            predicted_trend <- predict(trend_model, newdata = data.frame(valid_years = year_values))

            # Calculate the mean temperature for this day-hour combination (from valid data)
            mean_temp <- mean(valid_temp[valid_years %in% baseline_years], na.rm = TRUE)

            # Detrend by removing trend and adding back the mean
            values - predicted_trend + mean_temp
          }
        }
      }
    ) %>%
    ungroup()

  return(detrended_data$detrended)
}

# Apply the detrending function to your data
temp_detrended <- detrend_hourly_mc(climdata)
relhum_detrended <- detrend_hourly_mc(climdata, mc_col = "relhum")

# Add detrended columns to the original data frame
climdata_detrended <- climdata
climdata_detrended$relhum_detrended <- relhum_detrended
climdata_detrended$temp_detrended <- temp_detrended

# Summary statistics to verify the detrending
print("\nSummary of original vs detrended temperature:")
summary_stats <- climdata_detrended %>%
  summarise(
    original_temp_mean = mean(temp, na.rm = TRUE),
    original_temp_sd = sd(temp, na.rm = TRUE),
    detrended_temp_mean = mean(temp_detrended, na.rm = TRUE),
    detrended_temp_sd = sd(temp_detrended, na.rm = TRUE),
    mean_trend_removed = mean(temp_trend, na.rm = TRUE)
  )

print(summary_stats)

# Verify that long-term trend has been removed
# Calculate linear trend in annual means for both original and detrended data
trend_check <- annual_comparison %>%
  summarise(
    original_trend_slope = lm(original_annual_mean ~ year)$coefficients[2],
    detrended_trend_slope = lm(detrended_annual_mean ~ year)$coefficients[2]
  )

print("\nTrend verification (slope of annual means):")
print(paste("Original data trend:", round(trend_check$original_trend_slope, 4), "°C/year"))
print(paste("Detrended data trend:", round(trend_check$detrended_trend_slope, 4), "°C/year"))

# -- Plot data

climdata_detrended <- read_csv(file.path(out_dir, "climdata_era5_cmip6_1981-2100_ssp245_no_cc.csv"))

# ---- Annual statistics for temperature ----
annual_mc <- climdata_detrended %>%
  mutate(year = lubridate::year(obs_time)) %>%
  group_by(year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = TRUE),
    sd_temp   = sd(temp, na.rm = TRUE),
    mean_temp_detrended = mean(temp_detrended, na.rm = TRUE),
    sd_temp_detrended   = sd(temp_detrended, na.rm = TRUE),
    mean_relhum = mean(relhum, na.rm = TRUE),
    sd_relhum   = sd(relhum, na.rm = TRUE),
    mean_relhum_detrended = mean(relhum_detrended, na.rm = TRUE),
    sd_relhum_detrended   = sd(relhum_detrended, na.rm = TRUE)
  )

# ---- Plot for temperature ----
p_temp <- ggplot(annual_mc, aes(x = year)) +

  # Climate change
  geom_ribbon(
    aes(ymin = mean_temp - sd_temp,
        ymax = mean_temp + sd_temp,
        fill = "Climate change"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_temp, color = "Climate change"), size = 1) +
  geom_point(aes(y = mean_temp, color = "Climate change"), size = 1) +

  # Baseline (detrended)
  geom_ribbon(
    aes(ymin = mean_temp_detrended - sd_temp_detrended,
        ymax = mean_temp_detrended + sd_temp_detrended,
        fill = "Baseline"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_temp_detrended, color = "Baseline"), size = 1) +
  geom_point(aes(y = mean_temp_detrended, color = "Baseline"), size = 1) +

  scale_color_manual(
    values = c("Climate change" = "#F7766E",
               "Baseline" = "#004AAD")
  ) +
  scale_fill_manual(
    values = c("Climate change" = "#F7766E",
               "Baseline" = "#004AAD")
  ) +

  labs(
    x = "Year",
    y = "Temperature (°C)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 15)
  )

# ---- Plot for relative humidity ----
p_hum <- ggplot(annual_mc, aes(x = year)) +

  geom_ribbon(
    aes(ymin = mean_relhum - sd_relhum,
        ymax = mean_relhum + sd_relhum,
        fill = "Climate change"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_relhum, color = "Climate change"), size = 1) +
  geom_point(aes(y = mean_relhum, color = "Climate change"), size = 1) +

  geom_ribbon(
    aes(ymin = mean_relhum_detrended - sd_relhum_detrended,
        ymax = mean_relhum_detrended + sd_relhum_detrended,
        fill = "Baseline"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_relhum_detrended, color = "Baseline"), size = 1) +
  geom_point(aes(y = mean_relhum_detrended, color = "Baseline"), size = 1) +

  scale_color_manual(
    values = c("Climate change" = "#F7766E",
               "Baseline" = "#004AAD")
  ) +
  scale_fill_manual(
    values = c("Climate change" = "#F7766E",
               "Baseline" = "#004AAD")
  ) +

  labs(
    x = "Year",
    y = "Relative Humidity (%)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 15)
  )

# ---- Save to PDF ----
pdf(file.path(figs_dir, "compare_temp_relhum_cmip6_annual_cc_vs_no_cc_1981-2100_119ts_regua.pdf"),
    width = 10, height = 5)

(p_temp | p_hum) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()

# -- Save data
climdata_detrended["temp"] <- climdata_detrended["temp_detrended"]
climdata_detrended["relhum"] <- climdata_detrended["relhum_detrended"]
climdata_detrended <- climdata_detrended %>%
  select(-temp_detrended, -relhum_detrended)
write_csv(climdata_detrended, file.path(out_dir, "climdata_era5_cmip6_1981-2100_ssp245_no_cc.csv"))
