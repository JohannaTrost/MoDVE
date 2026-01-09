# Create species matrices
source("utils.R")

library("dplyr")

sampleEnvironment <- function(dimPlot, timeVals, numNiches) {
    
        # Define the ranges
        xVals <- 1:dimPlot[1]
        yVals <- 1:dimPlot[2]
        zVals <- 1:dimPlot[3]

        # Total number of possible combinations
        totalCombs <- length(xVals) * length(yVals) * length(zVals) * length(timeVals)

        # Sample unique indices
        sampleIndices <- sample(totalCombs, numNiches)

        # Convert sampled indices to actual combinations using array indexing
        getCombination <- function(index) {
          timeIndex <- ((index - 1) %% length(timeVals)) + 1
          zIndex <- (((index - 1) %/% length(timeVals)) %% length(zVals)) + 1
          yIndex <- ((((index - 1) %/% length(timeVals)) %/% length(zVals)) %% length(yVals)) + 1
          xIndex <- ((((index - 1) %/% length(timeVals)) %/% length(zVals)) %/% length(yVals)) + 1

          c(x = xVals[xIndex],
            y = yVals[yIndex],
            z = zVals[zIndex],
            time = timeVals[timeIndex])
        }

        # Apply the index-to-combination conversion
        result <- t(sapply(sampleIndices, getCombination))

        return(as.data.frame(result))
}


AgeMaturityMetabolic <- function(InterceptAgeMaturity, ScalingAgeMaturity, Mass) {
    return(InterceptAgeMaturity * (Mass^ScalingAgeMaturity))
}


Height2Light <- function(Height, kL, LAI, Imax) {
    # Convert height to light
    return(Imax * exp(-kL * LAI * (1 - Height)))
}


CheckNicheValues <- function(Niches) {
    # Humidity
    idx <- Niches$MinHum >= Niches$OptHum
    if (any(idx)) message(sum(idx), " case(s): MinHum >= OptHum, adjusting MinHum.")
    Niches[idx, "MinHum"] <- Niches[idx, "OptHum"] - 1

    idx <- Niches$MaxHum <= Niches$OptHum
    if (any(idx)) message(sum(idx), " case(s): MaxHum <= OptHum, adjusting MaxHum.")
    Niches[idx, "MaxHum"] <- Niches[idx, "OptHum"] + 1

    # Temperature
    idx <- Niches$MinTemp >= Niches$OptTemp
    if (any(idx)) message(sum(idx), " case(s): MinTemp >= OptTemp, adjusting MinTemp.")
    Niches[idx, "MinTemp"] <- Niches[idx, "OptTemp"] - 1

    idx <- Niches$MaxTemp <= Niches$OptTemp
    if (any(idx)) message(sum(idx), " case(s): MaxTemp <= OptTemp, adjusting MaxTemp.")
    Niches[idx, "MaxTemp"] <- Niches[idx, "OptTemp"] + 1

    # Wind
    idx <- Niches$MinWind >= Niches$OptWind
    if (any(idx)) message(sum(idx), " case(s): MinWind >= OptWind, adjusting MinWind.")
    Niches[idx, "MinWind"] <- Niches[idx, "OptWind"] - 1

    idx <- Niches$MaxWind <= Niches$OptWind
    if (any(idx)) message(sum(idx), " case(s): MaxWind <= OptWind, adjusting MaxWind.")
    Niches[idx, "MaxWind"] <- Niches[idx, "OptWind"] + 1

    # Light
    idx <- Niches$MinLight >= Niches$OptLight
    if (any(idx)) message(sum(idx), " case(s): MinLight >= OptLight, adjusting MinLight.")
    Niches[idx, "MinLight"] <- Niches[idx, "OptLight"] - 1

    idx <- Niches$MaxLight <= Niches$OptLight
    if (any(idx)) message(sum(idx), " case(s): MaxLight <= OptLight, adjusting MaxLight.")
    Niches[idx, "MaxLight"] <- Niches[idx, "OptLight"] + 1

    return(Niches)
}


# Parse input configuration file
config <- parse_config()

# ============================================================================
# RNG seed
seed <- config$seed
set.seed(seed, kind="Mersenne-Twister")  # integer for fixed seed or NULL for random

# Parameters that need to be specified/checked before running this script
MainOutputDirectory <- config$MainOutputDirectory
PathUniqueEnvVarCombs <- config$PathUniqueEnvVarCombs
Directorymicrohabitat <- config$Directorymicrohabitat

# Folder to save species trait matrices
# NameSpeciesPool <- "IntAgeMat_2_IntRec_70"  # Give meaningful name (the species type is automatically added to the name)

# Define number of species in species pool and total number of species pools to be created
numSpeciesPools <- config$numSpeciesPools
NumberOfSpecies <- config$NumberOfSpecies

# Time steps
initialTimeStep <- config$initialTimeStep
timeSteps <- config$timeSteps

# The following option defines if correlations between traits are consider or not
# CorrelationMassAgeOfMaturity <- 1  # Correlation between the mass and the age of maturity (this also influences the growth rate)
CorrelationMassRecruitment <- config$CorrelationMassRecruitment  # Correlation between the mass and the recruitment

InterceptAgeMaturity <- config$InterceptAgeMaturity
ScalingAgeMaturity <- config$ScalingAgeMaturity  # Scaling factor according to metabolic theory

# If correlations are choosen, the following parameters define the shape of the correlations
# 1. Correlations if CorrelationMassAgeOfMaturity=1
MaxMassRangeCorr <- config$MaxMassRangeCorr  # maximum mass of species/functional types (g)
# AgeAtMaturityRangeCorr=[1 15]  # age at which maturity is reaches (years) #Comment2019 => should not ne needed, delete
AgeAtMaturityDevCorr <- config$AgeAtMaturityDevCorr  # relative deviation from mean age of maturity

# 2. Correlations if CorrelationMassRecruitment=1
RecruitmentNormalizeAtSize1Corr <- config$RecruitmentNormalizeAtSize1Corr  # Factor converting the reproductive biomass to potential recruits
SlopeRecruitmentCorr <- 0 * RecruitmentNormalizeAtSize1Corr  # Slope of the correlation between mass and recruitment
RecruitmentInvestmentRelMeanCorr <- config$RecruitmentInvestmentRelMeanCorr  # Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
RecruitmentInvestmentRelDevCorr <- config$RecruitmentInvestmentRelDevCorr  # The relative deviation from the mean recruitment
RecruitmentIncMaxCorr <- config$RecruitmentIncMaxCorr

# Parameters of light model (needed to convert the height nicht and the light niche, these values do not have to
# be the same as used in the microhabitat matrices)
kL <- config$kL  # light extinction coefficient
Imax <- config$Imax  # maximum light intensity
LAI <- config$LAI  # leaf area index

# Define a species pool type
# SpeciesPoolType <- 0

# ============================================================================
# Define trait (ranges) if random species pool(SpeciesPoolType=0) is choose
# If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0), traits are randomly choosen from the following ranges
MaxMassRandom <- config$MaxMassRandom  # maximum mass of species/functional types (g)
MaxMassLogScaleRandom <- config$MaxMassLogScaleRandom  # define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturityRandom <- config$AgeAtMaturityRandom  # age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Random <- config$RecruitmentNormalizeAtSize1Random  # This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanRandom <- config$RecruitmentInvestmentRelMeanRandom  # Not that the effective recruitment is RecruitmentNormalizeAtSize1Random*RecruitmentInvestmentRelMeanRandom
RecruitmentIncRandom <- config$RecruitmentIncRandom
MassAtMaturityRelativeRandom <- config$MassAtMaturityRelativeRandom  # Relative mass in relation to maximum Size
DispersalKernelRandom <- config$DispersalKernelRandom  # The higher this values, the more local is the dispersal
DispersalKernelAsymmetryRandom <- config$DispersalKernelAsymmetryRandom  # The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)

# microclimate parameters
Random <- config$Random  # Flag to indicate if random MC niches
HeightBreadthRandom <- config$HeightBreadthRandom  # Relative height
LightBreadthRandom <- config$LightBreadthRandom
HumBreadthRandom <- config$HumBreadthRandom  # Relative humidity
TempBreadthRandom <- config$TempBreadthRandom  # Temperature in degrees Celsius
WindBreadthRandom <- config$WindBreadthRandom  # Wind speed in m/s

# If random light breadth is not defined, infer it from the random height breadth
if (length(LightBreadthRandom) == 0) {
    MinLightRandom <- Height2Light(HeightBreadthRandom[1], kL, LAI, Imax)
    MaxLightRandom <- Height2Light(HeightBreadthRandom[2], kL, LAI, Imax)
    LightBreadthRandom <- c(MinLightRandom, MaxLightRandom)
}

dimPlot <- readRDS(file.path(Directorymicrohabitat, "dimPlot.rds")) # Load plot dimensions

# Get microhabitat matrix variables - dynamic handling of selected variables
microhabitatVariableFlags <- config$microhabitatVariableFlags
microhabitat_var_names <- c(
  "TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
  "HumNicheOpt", "TempNicheOpt", "WindNicheOpt"
)
# Only keep active options
active_options <- microhabitat_var_names[microhabitatVariableFlags == 1]
# Assign indices
microhabitat_index_list <- setNames(seq_along(active_options), active_options)
# Check which env. variables are available and make named list with flags
env_var_names <- c("LightNicheOpt", "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
env_var_flags <- env_var_names %in% active_options
names(env_var_flags) <- c("Light", "Hum", "Temp", "Wind")

# where to save
SaveDirectory <- file.path(MainOutputDirectory)
dir.create(SaveDirectory, recursive=TRUE)

# Species trait matrix headers
ColumnHeaders <- c("SpeciesID", "MaximumMass", "MassAtMaturity", "GrowthRate",
               "DispersalKernel", "DispersalKernelAsymmetry", "RecruitmentInvestmentRel",
               "RecruitmentInc",
               "MaxRecruitsAtMaxMass", "MaxRecruitsAtMassAtMaturity", "AgeAtMaturity",
               "MinLight", "MaxLight", "OptimumLight",
               "LightResponseA", "LightResponseB", "LightResponseC",
               "MinHum", "MaxHum", "OptimumHum",
               "MinTemp", "MaxTemp", "OptimumTemp",
               "MinWind", "MaxWind", "OptimumWind",
               "DispersalKernelWindEffect"
)

# ===== Init random MC niches parameters =====

numNiches <- numSpeciesPools * NumberOfSpecies

# Store the breadth variables in a list
Breadths <- list(
  Hum = HumBreadthRandom,
  Temp = TempBreadthRandom,
  Wind = WindBreadthRandom,
  Light = LightBreadthRandom  # Convert relative light breadth to absolute values
)

# Generate random optimal environmental values
WindMargin <- runif(numNiches, 0.001, 1)
OptEnvVals <- data.frame(
  OptLight = runif(numNiches, min = Breadths$Light[1] + 1, max = Breadths$Light[2] - 1),
  OptHum = runif(numNiches, min = Breadths$Hum[1] + 1, max = Breadths$Hum[2] - 1),
  OptTemp = runif(numNiches, min = Breadths$Temp[1] + 1, max = Breadths$Temp[2] - 1),
  OptWind = runif(numNiches, min = Breadths$Wind[1] + WindMargin, max = Breadths$Wind[2] - WindMargin)
)

# Draw light breadth
possibleBreadths <- cbind(as.vector(OptEnvVals$OptLight) - Breadths$Light[1],
                          Breadths$Light[2] + 10 - as.vector(OptEnvVals$OptLight))
maxLightBreadths <- apply(possibleBreadths, 1, min)  # Ensure that light breadth does not exceed the range of light values
LightBreadths <- runif(numNiches, min = 1, max = maxLightBreadths)

Niches <- OptEnvVals %>%
  mutate(
    MinHum = runif(numNiches, min = Breadths$Hum[1], max = as.vector(OptHum)),
    MaxHum = runif(numNiches, min = as.vector(OptHum), max = Breadths$Hum[2]),
    MinTemp = runif(numNiches, min = Breadths$Temp[1], max = as.vector(OptTemp)),
    MaxTemp = runif(numNiches, min = as.vector(OptTemp), max = Breadths$Temp[2]),
    MinWind = runif(numNiches, min = Breadths$Wind[1], max = as.vector(OptWind)),
    MaxWind = runif(numNiches, min = as.vector(OptWind), max = Breadths$Wind[2]),
    MinLight = as.vector(OptLight) - LightBreadths,
    MaxLight = as.vector(OptLight) + LightBreadths
  )

# Ensure that min < opt < max -> very unlikely that min == opt or max == opt
Niches <- CheckNicheValues(Niches)

# Add species pool number and species number
Niches["SpeciesPool"] <- rep(seq_len(numSpeciesPools), each = NumberOfSpecies)
Niches["SpeciesID"] <- rep(seq_len(NumberOfSpecies), times = numSpeciesPools)

# Main loop (for random generation of species pool)
for (Num in seq_len(numSpeciesPools)) {

    # Trait matrix where the trait information of each species is saved
    SpeciesTraitMatrix <- matrix(0, NumberOfSpecies, length(ColumnHeaders))

    for (NumSpecies in seq_len(NumberOfSpecies)) {

        # ============================================================================
        # Maximum Size, size at maturity and growth rate (here, we are
        # choosing from log scale because usually, there are more smaller species than larger ones)
        if (MaxMassLogScaleRandom == 1) {
            MaxMassLog <- runif(1, min=log10(MaxMassRandom[1]), max=log10(MaxMassRandom[2]))
            MaxMass <- 10^MaxMassLog
        } else if (MaxMassLogScaleRandom == 0) {
            MaxMass <- runif(1, min=MaxMassRandom[1], max=MaxMassRandom[2])
        }

        # The age at maturity is the given by
        AgeAtMaturity <- AgeMaturityMetabolic(InterceptAgeMaturity, ScalingAgeMaturity, MaxMass)
        AgeAtMaturity <- AgeAtMaturity * runif(1, min=1-AgeAtMaturityDevCorr, max=1+AgeAtMaturityDevCorr)  # Add stochasticity

        # We assume that the mass at maturity is a function of the maximum size
        MassAtMaturity <- MaxMass * runif(1, min=MassAtMaturityRelativeRandom[1], max=MassAtMaturityRelativeRandom[2])

        # In the model, we are approximating growth by a Betalanffy growth curve, which generall is as follows
        # SizeFunctionOfAge=@(MaxMass,K,Age) (MaxMass*(1-exp(-K*(Age))));
        # By assuming that the Betalanffy growth curve crosses the point AgeAtMaturity/Size MaturityMassAtMaturity,
        # the growth rate K of this function can be calculated:
        K <- -(log(1) + log(1 - (MassAtMaturity / MaxMass))) / AgeAtMaturity

        # ============================================================================
        # Recruitment traits

        # 2. Correlations if CorrelationMassRecruitment=1
        # RecruitmentNormalizeAtSize1Corr <- 70  # Factor converting the reproductive biomass to potential recruits
        # SlopeRecruitmentCorr <- 0 * RecruitmentNormalizeAtSize1Corr  # Slope of the correlation between mass and recruitment
        # RecruitmentInvestmentRelMeanCorr <- 0.1  # Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
        # RecruitmentInvestmentRelDevCorr <- 0.25  # The relative deviation from the mean recruitment
        # RecruitmentIncMaxCorr <- 0

        if (CorrelationMassRecruitment == 1) {
            RecruitmentInvestmentRel <- runif(1, min=RecruitmentInvestmentRelMeanCorr * (1 - RecruitmentInvestmentRelDevCorr), max=RecruitmentInvestmentRelMeanCorr * (1 + RecruitmentInvestmentRelDevCorr))
            RecruitmentNormalizeAtSize1 <- RecruitmentNormalizeAtSize1Corr  # Factor converting the reproductive biomass to potential recruits
            SlopeRecruitment <- 0  # Slope of the correlation between mass and recruitment
            InterceptRecruitment <- RecruitmentNormalizeAtSize1
            RecruitmentInc <- 0
        } else if (CorrelationMassRecruitment == 0) {
            RecruitmentInvestmentRel <- runif(1, min=RecruitmentInvestmentRelMeanRandom[1], max=RecruitmentInvestmentRelMeanRandom[2])
            RecruitmentNormalizeAtSize1 <- runif(1, min=RecruitmentNormalizeAtSize1Random[1], max=RecruitmentNormalizeAtSize1Random[2])
            SlopeRecruitment <- 0  # No slope if no correlation is choosen
            InterceptRecruitment <- RecruitmentNormalizeAtSize1 - SlopeRecruitment
            RecruitmentInc <- runif(1, min=RecruitmentIncRandom[1], max=RecruitmentIncRandom[2])  # Not meaningful if no correlation
        }

        # Draw dispersal traits randomly
        DispersalKernel <- runif(1, min=DispersalKernelRandom[1], max=DispersalKernelRandom[2])
        DispersalKernelAsymmetry <- runif(1, min=DispersalKernelAsymmetryRandom[1], max=DispersalKernelAsymmetryRandom[2])
        DispersalKernelWindEffect <- runif(1, min = 0, max = 1)

        # ============================================================================
        # Traits of ecologcial light niche
        # Calculate parameters of parabolic response curve y=ax^2+bx+c
        # We assume that the function is a paraboloid which goes trough
        # three points (MinLight/0) (MaxLight/0) (OptimumLight/1)

        x1 <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MinLight"]
        x2 <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MaxLight"]
        x3 <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "OptLight"]

        y1 <- 0
        y2 <- 0
        y3 <- 1

        a <- (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / ((x1-x2) * (x1-x3) * (x3-x2))
        b <- (x1^2*(y2-y3) + x2^2*(y3-y1) + x3^2*(y1-y2)) / ((x1-x2) * (x1-x3) * (x2-x3))
        c <- (x1^2*(x2*y3 - x3*y2) + x1*(x3^2*y2 - x2^2*y3) + x2*x3*y1*(x2-x3)) / ((x1-x2) * (x1-x3) * (x2-x3))
        # In the model, based on these parameters the light response for each species can be calculated:
        # Parabol=@(a,b,c,x) a*x^2+b*x+c;


        # ============================================================================

        # Assign trait values for each species
        SpeciesTraitMatrix[NumSpecies, 1] <- NumSpecies  # Species or functional type
        SpeciesTraitMatrix[NumSpecies, 2] <- MaxMass  # Maximum mass
        SpeciesTraitMatrix[NumSpecies, 3] <- MassAtMaturity  # Mass at maturity
        SpeciesTraitMatrix[NumSpecies, 4] <- K  # Species-specfic growth rate
        SpeciesTraitMatrix[NumSpecies, 5] <- DispersalKernel  # Dispersal: factor b in negative exp funtion
        SpeciesTraitMatrix[NumSpecies, 6] <- DispersalKernelAsymmetry  # Dispersal: factor b in negative exp funtion
        SpeciesTraitMatrix[NumSpecies, 7] <- RecruitmentInvestmentRel  # anual reproductive allocation in relation to vegetative biomass
        SpeciesTraitMatrix[NumSpecies, 8] <- RecruitmentInc  # Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
        SpeciesTraitMatrix[NumSpecies, 9] <- (InterceptRecruitment)*RecruitmentInvestmentRel  # Potential maximum number of recruits at maximum mass
        SpeciesTraitMatrix[NumSpecies, 10] <- (InterceptRecruitment+SlopeRecruitment*MassAtMaturity)*RecruitmentInvestmentRel  # Potential maximum number of recruits at mass at maturity
        SpeciesTraitMatrix[NumSpecies, 11] <- AgeAtMaturity  # Average age at maturity under optimal conditions
        SpeciesTraitMatrix[NumSpecies, 12] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MinLight"]  # Min Light conditions
        SpeciesTraitMatrix[NumSpecies, 13] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MaxLight"]  # Max Light conditions
        SpeciesTraitMatrix[NumSpecies, 14] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "OptLight"]  # Optimum Light conditions
        SpeciesTraitMatrix[NumSpecies, 15] <- a  # Factor a of light response function
        SpeciesTraitMatrix[NumSpecies, 16] <- b  # Factor b of light response function
        SpeciesTraitMatrix[NumSpecies, 17] <- c  # Factor c of light response function
        SpeciesTraitMatrix[NumSpecies, 18] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MinHum"]  # Minimum relative humidity
        SpeciesTraitMatrix[NumSpecies, 19] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MaxHum"]  # Maximum relative humidity
        SpeciesTraitMatrix[NumSpecies, 20] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "OptHum"]  # Optimum relative humidity
        SpeciesTraitMatrix[NumSpecies, 21] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MinTemp"]  # Minimum temperature in degrees Celsius
        SpeciesTraitMatrix[NumSpecies, 22] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MaxTemp"]  # Maximum temperature in degrees Celsius
        SpeciesTraitMatrix[NumSpecies, 23] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "OptTemp"]  # Optimum temperature in degrees Celsius
        SpeciesTraitMatrix[NumSpecies, 24] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MinWind"]  # Minimum wind speed in m/s
        SpeciesTraitMatrix[NumSpecies, 25] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "MaxWind"]  # Maximum wind speed in m/s
        SpeciesTraitMatrix[NumSpecies, 26] <- Niches[Niches$SpeciesID == NumSpecies & Niches$SpeciesPool == Num, "OptWind"]  # Optimum wind speed in m/s
        SpeciesTraitMatrix[NumSpecies, 27] <- DispersalKernelWindEffect  # Wind effect on dispersal kernel [0, 1]
    }

    # ============================================================================
    # Save trait matrices

    # Create dataframe from matrix (including headers)
    SpeciesTraitMatrix_df <- as.data.frame(SpeciesTraitMatrix)
    names(SpeciesTraitMatrix_df) <- ColumnHeaders

    # Save trait dataframe
    SpeciesPoolFileName <- paste("SpeciesPool", Num, ".csv", sep="")
    write.csv(SpeciesTraitMatrix_df, file.path(SaveDirectory, SpeciesPoolFileName), row.names=FALSE)

    # ============================================================================
    # Save trait ranges used to generate the species pool
    TraitRanges <- matrix(0, 24, 2)

    TraitRanges[1, ] <- SlopeRecruitment
    TraitRanges[2, ] <- InterceptRecruitment
    TraitRanges[3, ] <- MaxMassRangeCorr
    # TraitRanges[4, ] <- AgeAtMaturityRangeCorr  # Comment2019 => should not ne needed, delete
    TraitRanges[5, ] <- AgeAtMaturityDevCorr
    TraitRanges[6, ] <- RecruitmentNormalizeAtSize1Corr
    TraitRanges[7, ] <- SlopeRecruitmentCorr
    TraitRanges[8, ] <- RecruitmentInvestmentRelMeanCorr
    TraitRanges[9, ] <- RecruitmentInvestmentRelDevCorr
    TraitRanges[10, ] <- RecruitmentIncMaxCorr
    TraitRanges[11, ] <- MaxMassRandom
    TraitRanges[12, ] <- AgeAtMaturityRandom
    TraitRanges[13, ] <- RecruitmentNormalizeAtSize1Random
    TraitRanges[14, ] <- RecruitmentInvestmentRelMeanRandom
    TraitRanges[15, ] <- RecruitmentIncRandom
    TraitRanges[16, ] <- MassAtMaturityRelativeRandom
    TraitRanges[17, ] <- DispersalKernelRandom
    TraitRanges[18, ] <- DispersalKernelAsymmetryRandom
    TraitRanges[19, ] <- MaxMassLogScaleRandom
    TraitRanges[20, ] <- HumBreadthRandom
    TraitRanges[21, ] <- TempBreadthRandom
    TraitRanges[22, ] <- WindBreadthRandom
    TraitRanges[23, ] <- LightBreadthRandom

    write.table(
        TraitRanges,
        file=file.path(SaveDirectory, "TraitRanges.csv"),
        sep=",",
        row.names=FALSE,
        col.names=FALSE,
    )

}
