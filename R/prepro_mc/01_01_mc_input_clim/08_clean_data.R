library(readr)

# --- Replace implausible values in climate data

# - No CC
climdata_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245_no_cc.csv"
climdata_reg <- read_csv(climdata_path)

climdata_reg[climdata_reg$relhum < 30, "relhum"] <- 30
climdata_reg[climdata_reg$relhum > 100, "relhum"] <- 100
climdata_reg[climdata_reg$temp > 36, "temp"] <- 36
climdata_reg[climdata_reg$temp < 8, "temp"] <- 8

new_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245_no_cc_v1.csv"
write_csv(climdata_reg, new_path)

# - With CC
climdata_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245.csv"
climdata_reg <- read_csv(climdata_path)

climdata_reg[climdata_reg$relhum < 30, "relhum"] <- 30
climdata_reg[climdata_reg$relhum > 100, "relhum"] <- 100
climdata_reg[climdata_reg$temp > 38, "temp"] <- 38
climdata_reg[climdata_reg$temp < 8, "temp"] <- 8

new_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245_v1.csv"
write_csv(climdata_reg, new_path)

check <- read_csv(new_path)
new_path <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua/scenarios/climdata_era5_cmip6_1981-2100_ssp245_no_cc_v1.csv"
check <- read_csv(new_path)