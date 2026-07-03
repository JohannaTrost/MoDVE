# ------
# Clean macroclimate data replacing implausible values

library(readr)

data_dir <- file.path("../modve_data")

# --- Replace implausible values in climate data

# - No CC
climdata_path <- file.path(data_dir, "mc_input", "regua", "scenarios", "climdata_era5_cmip6_1981-2100_ssp245_no_cc.csv")
climdata_reg <- read_csv(climdata_path)

climdata_reg[climdata_reg$relhum < 30, "relhum"] <- 30
climdata_reg[climdata_reg$relhum > 100, "relhum"] <- 100
climdata_reg[climdata_reg$temp > 50, "temp"] <- 50
climdata_reg[climdata_reg$temp < 8, "temp"] <- 8

new_path <- file.path(data_dir, "mc_input", "regua", "scenarios", "climdata_era5_cmip6_1981-2100_ssp245_no_cc_v1.csv")
write_csv(climdata_reg, new_path)

# - With CC
climdata_path <- file.path(data_dir, "mc_input", "regua", "scenarios", "climdata_era5_cmip6_1981-2100_ssp245.csv")
climdata_reg <- read_csv(climdata_path)

climdata_reg[climdata_reg$relhum < 30, "relhum"] <- 30
climdata_reg[climdata_reg$relhum > 100, "relhum"] <- 100
climdata_reg[climdata_reg$temp > 50, "temp"] <- 50
climdata_reg[climdata_reg$temp < 8, "temp"] <- 8

new_path <- file.path(data_dir, "mc_input", "regua", "scenarios", "climdata_era5_cmip6_1981-2100_ssp245_v1.csv")
write_csv(climdata_reg, new_path)
