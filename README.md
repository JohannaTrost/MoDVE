# MoDVE: An integrated approach to modeling microclimate niches of vascular epiphytes
A simulation pipeline for agent-based mechanistic modeling of epiphyte demography in a neotropical lowland forest. The repository includes a step by step pipeline, for simulating microclimate with Micropoint (Maclean, 2025), integrating forest structure and light conditions from forest model outputs of MoF3D (Petter et al., 2021a). It builds upon the epiphyte model by Petter et al. (2021b) introducing temperature and humidity niche axes.

Key components: 
* Modeling pipeline (MoDVE, Micropoint)
* Analyses 

## Overview
This project aims at providing insights into vascular epiphyte community dynamics in response to changes in relative humidity and temperature. Comparing a non-climate-change with the SSP2-4.5 scenario sheds light on a potential future pathway for canopy-dwelling species in terms of species richness and species distribution along the vertical forest profile, which can provide relatively steep relative humidity and temperature gradients. 

This repository provides a full model pipeline and the following simulation analyses and experiments: 
* Degenerate sensitivity analysis
* Comparison to the mid-domain effect expectation of species richness
* Comparison between a non-climate-change and climate change scenario
* Analysis of interspecific competition effects along the forest profile.

## Table of Contents
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [License](#license)

## Installation
For running the modeling pipline an R version >= 4.3.1 is required. 
For downloading CMIP6 data and for several post-simulation scripts python >= 3.11 is required. 

1. Clone the repository: git clone https://github.com/JohannaTrost/MoDVE.git
2. Then navigate to the repository folder: cd MoDVE
3. Open a R console and install required R libraries:
```r   
install.packages("devtools")
devtools::install_deps(".")
```
4. Open a terminal and install Python packages:
```bash
pip install . -r requirements.txt
```

## Project Structure

```bash
MoDVE/
├── DESCRIPTION
├── README.md
├── analysis # Simulation analyses
│   ├── climate_change
│   └── sensitivity
├── model_pipeline # Core simulation pipeline (run sequentially)
│   ├── 01_generate_microhabitat.R
│   ├── 02_simulate_microclimate
│   ├── 03_add_microclimate_dimensions.R
│   ├── 04_create_species_pools.R
│   ├── 05_compute_env_suitability.R
│   ├── 06_create_initial_distributions.R
│   ├── 07_run_model.R
│   ├── 08_replicate_diagnostic.R # Evaluate the no. replicates
│   └── utils.R
└── requirements.txt
```

The `model_pipeline/` scripts are designed to be executed sequentially as a command-line workflow (see [Usage](#usage)). Shared utility functions used across the pipeline are provided in `model_pipeline/utils.R`.

The `analysis/` folder contains multiple simulation analysis scripts requiring outputs from `model_pipeline/07_run_model.R`. They are divided into two branches:

#### Sensitivity analysis (`analysis/sensitivity/`)

A degenerate sensitivity analysis examining how varying microclimate gradients affect simulated diversity and vertical distribution patterns. Scripts proceed from varying input gradients (`01_vary_mc_gradients.R`), computing richness and species distribution metrics (`02_compute_richness.py`, `03_compute_species_distribution.R`) from model outputs, to calculating relative changes across scenarios (`04_relative_changes.R`) and visualizing differences (`05_plot_scenario_differences.py`).

#### Climate change analysis (`analysis/climate_change/`)

Compares the non-climate-change baseline against the SSP2-4.5 scenario across five modules:

| Folder | Content |
|---|---|
| `01_community/` | Community-level visualizations: demographics, vertical distribution shifts, abundance, and species richness |
| `02_mixed_effect_models/` | Mixed-effect models for position, range, and vertical richness and abundance responses |
| `03_traits_and_niches/` | Functional trait and niche analyses including changes in realized and potential niches, environmental limits, and niche filling figures |
| `04_mid_domain_effect/` | Comparison of simulated richness patterns against mid-domain effect expectations |
| `05_experiment_upward_shift/` | Mechanistic experiment isolating upward-shifting species: filtering, distribution setup, and visualization of resulting position shifts |

Note that scripts in folders `02` to `05` may depend on files generated in `01`. 

## Usage

### Model pipeline walkthrough

Here, we demonstrate how to use the model pipeline using the example of simulating epiphytes over 20 years at the forest site REGUA using the SSP2-4.5 climate scenario based on CMIP6 data.
To reproduce all simulations conducted for the publication, please use the configuration files provided in `modve_data_zenodo/cfgs`. 

#### 0. Download input data from Zenodo
First, download the necessary input data from Zenodo: 
```bash
TODO
```

#### 1. Generate microhabitat
Then, generate the microhabitat matrices from the MoF3D forest output or find the corresponding files in `../modve_data_zenodo/modve_output/regua/microhabitat`.
Open a terminal and `cd` into the project folder (MoDVE). Then, execute the script:
```bash
Rscript model_pipeline/01_generate_microhabitat.R --config ../modve_data_zenodo/cfgs/01_config.toml
```
This produces 4D matrices (xDim X yDim X zDim X nVariables) for each time step in `../modve_output/regua/microhabitat`: `microhabitatMatrix80.rds`, ..., `microhabitatMatrix100.rds`, and forest parameter files (`Forest_param_global.txt`, `Forest_param_pass0.txt`, `dimPlot.rds`).

#### 2. Simulate microclimate

In the following scripts, any configurable variables and parameters can be found preceeding the initial library imports. This step can also be skipped, as the resulting microclimate simulations are also available on Zenodo (TODO link) in the directory `modve_data_zenodo/mc_output`.

##### 2.1 Climate input

1. <ins>Download CMIP6 data:</ins>
   * Register at the CEDA Archive ![https://catalogue.ceda.ac.uk/uuid/c107618f1db34801bb88a1e927b82317/](https://catalogue.ceda.ac.uk/uuid/c107618f1db34801bb88a1e927b82317/)
   * Generate an access token
   * Execute CMIP6 downloader for historical and SSP2-4.5 data:
  ```bash
cd model_pipeline/02_simulate_microclimate/01_climate_inputs
OUTDIR="../../../../modve_data/mc_input/climate/cmip6_ceda"

python 01_cmip6_downloader.py --help
python 01_cmip6_downloader.py -f 1981 -l 2014 -o $OUTDIR --scenario historical --token-file <path/to/token.txt>
python 01_cmip6_downloader.py -f 2015 -l 2100 -o $OUTDIR --scenario ssp245 --token <your token>
  ```
This will produce the files `baf_ensemble_day_<scenario>_<year>.nc` and `baf_<var>_<scenario>_<year>.nc` with `year` 1981-2100, `scenario` either historical or ssp245, and `var` including daily precipitation (pr_day), wind (sfcWind), relative humidity (hurs), pressure (ps_day), and temperature (tas).

2. <ins>Download ERA5 data:</ins>
```bash
CLIMDIR="../../../../modve_data/mc_input/climate"
for year in {1981..2025}; do
  Rscript 02_era5_downloader.R $CLIMDIR/era5_raw year
end
```
3. <ins>Process ERA5 data:</ins>
```bash
Rscript 03_prepro_era5.R $CLIMDIR/era5_raw $CLIMDIR/era5_processed $(seq 1981 2025)
```
4. <ins>Execute R scripts</ins> step by step to select a subregion and merge ERA5 and CMIP6 data for complete hourly variables (`04_subset_region_and_merge.R`), detrend the climate change data (SSP2-4.5) for a baseline scenario (`05_generate_baseline_climate.R`) and clean the resulting data sets (`06_clean_data.R`). Note that any configurable variables will be listed after loading required libraries.
5. Go back to the project folder directory:
   ```bash
   cd ../../../
   ```

##### 2.2 Vegetatation, soil, and albedo

1. <ins>Vegetation:</ins>
  * Generate a Google Earth Enginge (GEE) project, initialize it and configure the respective variables in the script `model_pipeline/02_simulate_microclimate/02_surface_inputs/01_vegetation.R`
  * Walk through the script downloading landcover and vegetation height with GEE and generating vegeation inputs.
2. <ins>Albedo:</ins> Use the script `model_pipeline/02_simulate_microclimate/02_surface_inputs/02_albedo.R` to download and visualize MODIS RGB and CIR data stored in `../modve_data/mc_input/albedo` and `../modve_figs`.
3. <ins>Soil:</ins> Generate soil and leaf reflectance data and rasters with soil parameters with `model_pipeline/02_simulate_microclimate/02_surface_inputs/03_soil.R`
4. Crop soil and vegetation to 50 m x 50 m rasters for either the REGUA or Pirineus forest stand in the Atlatnic forest and process list of rasters to specific format required by Micropoint with `model_pipeline/02_simulate_microclimate/02_surface_inputs/04_format_subregion.R`
5. Replace vegetation height and PAI with data from MoF3D simulations with `model_pipeline/02_simulate_microclimate/02_surface_inputs/05_mc_input_mof3d.R` 

##### 2.3 Simulation with Micropoint 

1. Simulate microclimate for each year (in this example 1981-2000, for SSP2-4.5) processing a chunk of 50 cells at a time (for a total of 2500 cells):
```bash
#!/bin/bash

# Change directory to simulation folder
cd model_pipeline/02_simulate_microclimate/03_simulation

# Path to simulation config file
cfg="../../../../modve_data_zenodo/cfgs/02_1_mc_sim.toml"

# Loop over cells in steps of 50 (2500 cells, 50 cells in parallel)
for chunk in $(seq 1 50 2501); do
  # Loop over years to simulate
  for year in {1981..2000}; do
    # Get time step corresponding to MoF3D simulation steps
    t=$((year - 1901))

    # Run the simulation script
    Rscript 01_simulate_mc.R --config $cfg --year "$year" --timestep "$t" --chunk "$chunk"
  done
done
```
2. Merge all cells into a microclimate matrix:
```bash
# Path to config file
cfg="../../../../modve_data_zenodo/cfgs/02_2_mc_process.toml"

# Loop over years to simulate
for year in {1981..2000}; do
  # Get time step corresponding to MoF3D simulation steps
  t=$((year - 1901))

  # Run the simulation script
  Rscript 02_merge_cells.R --config $cfg --year "$year" --timestep "$t"
done
```

##### 2.4 Visualize and validate microclimate simulations

The following scripts in `model_pipeline/02_simulate_microclimate/04_diagnostics` can be used to reproduce the comparison between simulated and measured microclimate including computing errors, correlations and plotting time series and vertical profiles. For this, make sure that the project data have been downloaded from Zenodo (see step 0.) and copied into the parent directory of the project root or adjust paths manually in the respective scripts. 

```bash
└── 04_diagnostics
    ├── validate_mc_empirical_2024.R
    ├── validate_mc_gradient_empirical_2025.R
    ├── visualize_mc_pirineus_1906-2024_era5.R
    ├── visualize_mc_regua_1981-2100_ssp245_no_cc.R
    └── visualize_vertical_mc_regua_ssp245_no_cc.rmd
```

#### 3. Add microclimate to microhabitat matrices
Add the simulated microclimate variables to the microhabitat matrices or find the corresponding files in `../modve_data_zenodo/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/microhabitat_mc/forest0`. 
```bash
Rscript model_pipeline/03_add_microclimate_dimensions.R --config ../modve_data_zenodo/cfgs/03_config.toml
```
New `MicrohabitatMatrix<t>.rds` files will be in `../modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/microhabitat_mc/forest0` (with `t` = time step) and include matrices for REGUA, forest replicate 0, under the SSP 2-4.5 scenario (for other replicates and scenarios copy and adjust the config file). 

#### 4. Draw species
Next, generate the species for each species pool. To reproduce our results please use the species pools in `../modve_data_zenodo/modve_output/regua/species_pools`.

```bash
Rscript model_pipeline/04_create_species_pools.R --config ../modve_data_zenodo/cfgs/04_config.toml
```

#### 5. Precompute environmental suitability
First, compute environmental suitability scores (0-1) for each scenario, species, voxel, forest and time step: 

```bash
Rscript model_pipeline/05_compute_env_suitability.R --config ../modve_data_zenodo/cfgs/05_config.toml
```

Suitability scores will be stored in `../modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/suitability_scores/forest0/EnvSuitability/ID_SpeciesP_<pool>_TimeStep<t>.h5` (with `t` = time step, `pool` = no. species pool).
Second, scale the suitability scores across space, time, scenario, and forest for each time step (using the additional `singleStep` argument):

```bash
for singleStep in {80..100}; do
  Rscript model_pipeline/05_compute_env_suitability.R --config ../modve_data_zenodo/cfgs/05_config.toml *singleStep*
end
```

This will first compute the global maximum suitabilities for each species and species pool and save it in `../modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/suitability_scores/forest0/EnvSuitability/GlobalMaxSuitability_<pool>.h5`. Then, the suitability scores will be scaled for the given time step (`singleStep`) and written to `../modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/suitability_scores/forest0/EnvSuitability/ScaledSuitability_<pool>TimeStep<t>.h5`. The corresponding file with unscaled scores will be deleted thereafter.

#### 6. Initialize the spatial epiphyte distribution
For the starting point of the simulation, distribute epiphyte individuals for all species across the 3D forest stand according to suitable environmental factors. To reproduce our results please use the distributions in `../modve_data_zenodo/distribution`.

```bash
Rscript model_pipeline/06_create_initial_distributions.R --config ../modve_data_zenodo/cfgs/06_config.toml
```

The initial distribution will be saved in `../modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/distribution/forest0/ID_SpeciesP_<pool>_Rep_<replicate>.csv` (with `pool` = no. species pool, `replicate` = no. replicate simulation of the species pool, i.e. the same species pool with different random initial distributions).

#### 7. Simulate epiphyte communities
Finally, simulate epiphyte communities over 20 time steps (1980 ti 2000):

```bash
Rscript model_pipeline/07_run_model.R --config ../modve_data_zenodo/cfgs/07_config.toml
```

Epiphyte communities will be saved in `../modve_output/communities` and include epiphyte matrices for each time step t (`IndividualMatrixTimeStep<t>.csv`), summary statistics for each species and time step (`SpeciesSummary.csv`), summary statistics for the entire community (`CommunitySummary.csv`) and optionally the random number generator state (`random_state_seed.RData`).

#### 8. Evaluate the number of replicates used 
Evaluate whether the number of replicates (forests and species pools) is sufficient using a bootstrapping approach:

```bash
Rscript model_pipeline/08_replicate_diagnostic.R
```

How it works:
1. Draw an increasing number of replicates (1 to 30)
2. Compute the coefficient of variation (CV) across richness and abundance of these replicates
3. Repeat 100 times
4. Plot the number of draws against the CV (saved in `../modve_figs/climdata_era5_cmip6_1981-2100_ssp245`)

## License
This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

