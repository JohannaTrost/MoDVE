# MoDVE

MoF3D generates the 3D forest for epiphytes to inhabit, growing across multiple timesteps.

MoDVE scripts:

A1- converts the MoF3D output into readable microhabitat matrices for each timestep for the epiphytes to inhabit. The epiphytes depend on the amount of substrate (i.e. branch) and the amount of light in each voxel of the microhabitat matrix.

A2- generates several pools of epiphyte species with random traits

A3- distributes these randomly generated epiphytes in forest in the initial timestep.

A4- simulates these communities of random species growing in the forest microhabitats in the subsequent timesteps

B1- selects for viable/realistic species to include in the final simulation- i.e. species that won't heavily dominate the community, or rapidly go extinct. It then generates new pools of these realistic species to use.

B2- distributes these selected species in the forest in the initial timestep.

B3- simulates these communities growing the forest microhabitats of the subsequent timesteps.

# R Scripts

Install the following packages:

```bash
Rscript -e 'install.packages("optparse", repos="http://cran.uk.r-project.org")'
Rscript -e 'install.packages("configr", repos="http://cran.uk.r-project.org")'
```

## A1

Create a file named `config.toml`. Use the following content as a template, replacing the placeholder values with your data:

```toml
# This parameter determines which type of microhabitat matrices are generated:
# 1: real GroIMP forest with dynamics
# 2: static GroIMP forest (only forest at timeStepStart is used)
MicrohabitatType = 1

# Parameters of light model
# Light extinction coefficient
kL = 0.6
# How many rings around focal voxel to consider in light model
# (5 voxels in x and y direction)
DistVoxToConsider = 8

# Choose the forest parameters that shall be calculated and stored in the microhabitat matrix
# (this list can be extended for possible new applications of the epiphyte model)
# 1: use this variable
# 0: do not use it
TotalSurfaceAreaOpt = 1
SurfaceAreaLossOpt = 1
LightConditionsOpt = 1
AverageWeightedAngles = 0

# Parameters that need to be specified if MicrohabitatType=1 or MicrohabitatType=2
# Directory of GroIMP files (this directory is stored in the Microhabitat folder so that the
# connection to the input GroIMP files is always clear)
DirectoryGroIMP = "path/to/GroIMP/output/dir"
# Directory to save results
DirectorySaveMain = "path/to/output/dir"

ReplicateForest = 0

# start and end timestep
timeStepStart = 1
timeStepEnd = 40
```

Run the following command, replacing `path/to/toml` with the actual path to your `config.toml` file:

```bash
Rscript A1.R -i "path/to/toml"
```

## A2

Create a file named `config.toml`. Use the following content as a template, replacing the placeholder values with your data:

```toml
# Seed for random number generator (integer, optional)
# Comment it out to use a random seed instead
#seed = 42

MainOutputDirectory = "path/to/output/dir"

# Define number of species in species pool and total number of species pools to be created
numSpeciesPools = 100
NumberOfSpecies = 100

# The following option defines if correlations between traits are consider or not
CorrelationMassRecruitment = 1  # Correlation between the mass and the recruitment

InterceptAgeMaturity = 2
ScalingAgeMaturity = 0.25  # Scaling factor according to metabolic theory

# If correlations are choosen, the following parameters define the shape of the correlations
# 1. Correlations if CorrelationMassAgeOfMaturity=1
MaxMassRangeCorr = [2, 3000]  # maximum mass of species/functional types (g)
AgeAtMaturityDevCorr = 0.25  # relative deviation from mean age of maturity

# 2. Correlations if CorrelationMassRecruitment=1
RecruitmentNormalizeAtSize1Corr = 70  # Factor converting the reproductive biomass to potential recruits
RecruitmentInvestmentRelMeanCorr = 0.1  # Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
RecruitmentInvestmentRelDevCorr = 0.25  # The relative deviation from the mean recruitment
RecruitmentIncMaxCorr = 0

# Parameters of light model (needed to convert the height nicht and the light niche, these values do not have to
# be the same as used in the microhabitat matrices)
kL = 0.6  # light extinction coefficient
Imax = 900  # maximum light intensity
LAI = 6  # leaf area index

# Define trait (ranges) if random species pool(SpeciesPoolType=0) is choose
# If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0), traits are randomly choosen from the following ranges
MaxMassRandom = [2, 3000]  # maximum mass of species/functional types (g)
MaxMassLogScaleRandom = 1  # define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturityRandom = [1, 1]  # age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Random = [1, 20]  # This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanRandom = [0.07, 0.12]  # Not that the effective recruitment is RecruitmentNormalizeAtSize1Random*RecruitmentInvestmentRelMeanRandom
RecruitmentIncRandom = [0, 0]
MassAtMaturityRelativeRandom = [0.5, 0.7]  # Relative mass in relation to maximum Size
HeightBreadthRandom = [0.15, 0.7]  # Relative height
DispersalKernelRandom = [0.03, 0.5]  # The higher this values, the more local is the dispersal
DispersalKernelAsymmetryRandom = [0.5, 0.95]  # The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)
```

Run the following command, replacing `path/to/toml` with the actual path to your `config.toml` file:

```bash
Rscript A2.R -i "path/to/toml"
```

## A3

Create a file named `config.toml`. Use the following content as a template, replacing the placeholder values with your data:

```toml
# Seed for random number generator (integer, optional)
# Comment it out to use a random seed instead
#seed = 42

# Name of epiphyte model
# FolderEpiphyteModel='EM_20160213'; %I should think about naming
SingleSpeciesModel = 0  # 1: Single species model, 0: Community model

# Name of species pool
FolderSpeciesPools = "SP_Random_IntAgeMat_2_IntRec_70_TraitCorrOn"

# Directory where model is save and directory where microhabitat matrices
# are stored
DirectoryModelMain = "path/to/output"
DirectoryMicrohabitatMain = "path/to/microhabitat"
DirectorySpeciesPoolsMain = "path/to/species"

Replicate = 0

# Choose species pools to use and number of replicates per species pool
numSpeciesPools = [99, 100]  # Start and end number of  species pools
replicatePerSpeciesPool = 1  # Number of replicates per species pool
TimeStep = 1  # Time step for which the Initial distribution is generated

# The suitable voxel can either be the voxel
# with the highest available surface area (MethodVoxel=1)
# or a random voxel (MethodVoxel=0)
MethodVoxel = 0

# Define how many individuals per species are used, and how many of them are initially mature
# This variable defines if the NumberSpecies are total numbers irrespective
# of the model area (ScalingPerHa=0), or if the NumberSpecies or given per
# hectar and are scaled to the model area (ScalingPerHa=1)
ScalingPerHa = 0
IndividualsPerSpecies = 100
PercentageMaturePerSpecies = 50

# This parameter set the scaling between the
SurfaceBiomassScaling = 100  # cm^2 per m^2
Imax = 900
```

Run the following command, replacing `path/to/toml` with the actual path to your `config.toml` file:

```bash
Rscript A3.R -i "path/to/toml"
```
