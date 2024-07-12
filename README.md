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
seed = 42
```

Run the following command, replacing `path/to/toml` with the actual path to your `config.toml` file:

```bash
Rscript A2.R -i "path/to/toml"
```
