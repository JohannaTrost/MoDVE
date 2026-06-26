library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)


region <- "pirineus"
in_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/mc_input", region)
out_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/mc_input", region, "scenarios")

for (year in 1981:2024) {
  file <- paste(in_dir, paste0("climdata_era5_cmip6_", year, "_v3.csv"), sep = "/")
  climdata_year <- read_csv(file)

  if (year == 1981) {
    climdata_all <- climdata_year
  } else {
    climdata_all <- bind_rows(climdata_all, climdata_year)
  }
}

climdata_all <- climdata_all %>% select(-windspeed_anom)

write_csv(climdata_all, file.path(out_dir, "climdata_era5_cmip6_1981-2024_ssp245.csv"))

# --- 1. Scenario: 119 years until 2024

n_years_missing <- 119 - (2024 - 1981 + 1)

# Repeat first 5 years 15 times adding 75 to the time series to get 119 years
climdata_present_ssp245 <- climdata_all %>%
  filter(year(obs_time) %in% 1981:1985) %>%
  mutate(
    # Store original row count before expansion
    original_rows = n()
  ) %>%
  slice(rep(1:n(), times = 16)) %>%
  mutate(
    # Create sequence of target years (1911-1985)
    rep_group = rep(1:16, each = first(original_rows)),
    target_year_base = 1906 + (rep_group - 1) * 5,
    original_year = year(obs_time),
    target_years = target_year_base + (original_year - 1981),

    # Replace the year component while keeping month, day, hour, etc.
    obs_time = make_datetime(
      year = target_years,
      month = month(obs_time),
      day = day(obs_time),
      hour = hour(obs_time),
      min = minute(obs_time),
      sec = second(obs_time)
    )
  ) %>%
  select(-original_rows, -rep_group, -target_year_base, -original_year, -target_years)

# test
test <- climdata_present_ssp245 %>%
  mutate(year = year(obs_time), month = month(obs_time)) %>%
  select(year, month, obs_time) %>%
  group_by(year, month) %>%
  slice_head(n = 1)   # keeps only the first row per year-month

# Now combine climdata_present_ssp245 with climdata_all (from 1986 to 2024)
climdata_present <-
  bind_rows(climdata_present_ssp245,
            climdata_all %>% filter(year(obs_time) >= 1986)) %>%
  arrange(obs_time)

# Save scenario
write_csv(climdata_present, file.path(out_dir, "climdata_era5_cmip6_1906-2024_ssp245_119ts_v1.csv"))

# --- PLot CHECK

# Aggregate hourly to annual mean temperature
annual_temp <- climdata_present %>%
  mutate(year = year(obs_time)) %>%
  group_by(year) %>%
  summarise(annual_temp = mean(temp, na.rm = TRUE))

# Time series plot
p <- ggplot(annual_temp, aes(x = year, y = annual_temp)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 1) +
  labs(title = "Annual Mean Temperature",
       x = "Year",
       y = "Temperature (°C)") +
  theme_minimal(base_size = 14)

# Display the plot
pdf("../../figs/mc_input/compare_temp_cmip6_annual_present_cc_119ts_pirineus_v1.pdf")
print(p)
dev.off()

# --- 2. Scenario: 119 years until 2100 with CC - REGUA

region <- "regua"
in_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/mc_input", region)
out_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/mc_input", region, "scenarios")

for (year in 1981:2100) {
  file <- paste(in_dir, paste0("climdata_era5_cmip6_", year, "_v3.csv"), sep = "/")
  climdata_year <- read_csv(file)

  if (year == 1981) {
    climdata_all <- climdata_year
  } else {
    climdata_all <- bind_rows(climdata_all, climdata_year)
  }
}

climdata_all <- climdata_all %>% select(-windspeed_anom)

# Check if there are NAs
sum(is.na(climdata_all))

# --- PLot CHECK

# ---- Annual statistics for temperature ----
annual_temp <- climdata_all %>%
  mutate(year = lubridate::year(obs_time)) %>%
  group_by(year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = TRUE),
    sd_temp   = sd(temp, na.rm = TRUE)
  )

# ---- Annual statistics for relative humidity ----

annual_hum <- climdata_all %>%
  mutate(year = lubridate::year(obs_time)) %>%
  group_by(year) %>%
  summarise(
    mean_hum = mean(relhum, na.rm = TRUE),
    sd_hum   = sd(relhum, na.rm = TRUE)
  )

# ---- Plot for temperature ----
p_temp <- ggplot(annual_temp, aes(x = year, y = mean_temp)) +
  geom_ribbon(aes(ymin = mean_temp - sd_temp,
                  ymax = mean_temp + sd_temp),
              fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 1) +
  labs(x = "Year",
       y = "Temperature (°C)") +
  theme_minimal(base_size = 14)

# ---- Plot for relative humidity ----
p_hum <- ggplot(annual_hum, aes(x = year, y = mean_hum)) +
  geom_ribbon(aes(ymin = mean_hum - sd_hum,
                  ymax = mean_hum + sd_hum),
              fill = "darkgreen", alpha = 0.2) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "darkgreen", size = 1) +
  labs(x = "Year",
       y = "Relative Humidity (%)") +
  theme_minimal(base_size = 14)

# ---- Combine plots with patchwork ----
combined_plot <- p_temp / p_hum  # stacked vertically

# ---- Save to PDF ----
pdf("../../figs/mc_input/compare_temp_hum_cmip6_annual_cc_1981-2100_119ts_regua_v1.pdf",
    width = 8, height = 8)
print(combined_plot)
dev.off()


write_csv(climdata_all, file.path(out_dir, "climdata_era5_cmip6_1981-2100_ssp245.csv"))
