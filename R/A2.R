# Create species matrices

# ============================================================================
# Parameters that need to be specified/checked before running this script
MainOutputDirectory <- "/PATH/TO/OUTPUT"

# Folder to save species trait matrices
NameSpeciesPool <- "IntAgeMat_2_IntRec_70"  # Give meaningful name (the species type is automatically added to the name)

# Define number of species in species pool and total number of species pools to be created
numSpeciesPools <- 100
NumberOfSpecies <- 100

# The following option defines if correlations between traits are consider or not
CorrelationMassAgeOfMaturity <- 1  # Correlation between the mass and the age of maturity (this also influences the growth rate)
CorrelationMassRecruitment <- 1  # Correlation between the mass and the recruitment

InterceptAgeMaturity <- 2
ScalingAgeMaturity <- 0.25  # Scaling factor according to metabolic theory

# If correlations are choosen, the following parameters define the shape of the correlations
# 1. Correlations if CorrelationMassAgeOfMaturity=1
MaxMassRangeCorr <- c(2, 3000)  # maximum mass of species/functional types (g)
# AgeAtMaturityRangeCorr=[1 15]  # age at which maturity is reaches (years) #Comment2019 => should not ne needed, delete
AgeAtMaturityDevCorr <- 0.25  # relative deviation from mean age of maturity

# 2. Correlations if CorrelationMassRecruitment=1
RecruitmentNormalizeAtSize1Corr <- 70  # Factor converting the reproductive biomass to potential recruits
SlopeRecruitmentCorr <- 0 * RecruitmentNormalizeAtSize1Corr  # Slope of the correlation between mass and recruitment
RecruitmentInvestmentRelMeanCorr <- 0.1  # Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
RecruitmentInvestmentRelDevCorr <- 0.25  # The relative deviation from the mean recruitment
RecruitmentIncMaxCorr <- 0

# Parameters of light model (needed to convert the height nicht and the light niche, these values do not have to
# be the same as used in the microhabitat matrices)
kL <- 0.6  # light extinction coefficient
Imax <- 900  # maximum light intensity
LAI <- 6  # leaf area index

# Define a species pool type
SpeciesPoolType <- 0

# ============================================================================
# Define trait (ranges) if random species pool(SpeciesPoolType=0) is choose
# If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0), traits are randomly choosen from the following ranges
MaxMassRandom <- c(2, 3000)  # maximum mass of species/functional types (g)
MaxMassLogScaleRandom <- 1  # define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturityRandom <- c(1, 1)  # age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Random <- c(1, 20)  # This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanRandom <- c(0.07, 0.12)  # Not that the effective recruitment is RecruitmentNormalizeAtSize1Random*RecruitmentInvestmentRelMeanRandom
RecruitmentIncRandom <- c(0, 0)
MassAtMaturityRelativeRandom <- c(0.5, 0.7)  # Relative mass in relation to maximum Size
HeightBreadthRandom <- c(0.15, 0.7)  # Relative height
DispersalKernelRandom <- c(0.03, 0.5)  # The higher this values, the more local is the dispersal
DispersalKernelAsymmetryRandom <- c(0.5, 0.95)  # The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)

# ============================================================================
# Define which trait is varied if a sequential species pool(SpeciesPoolType=1) is choose.
# For proper results, only one of the following traits should be defined as sequence, while for the other traits, invariable trait values should to be specified.
# If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0),the following traits are used
MaxMassSeq <- c(2, 3000)  # maximum mass of species/functional types (g)
MaxMassLogScaleSeq <- 1  # define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturitySeq <- 2  # if no correlation is defined, this value is used
RecruitmentNormalizeAtSize1Seq <- 20  # This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanSeq <- 0.1  # anual investment in reproduction in relation to vegetative biomass
RecruitmentIncSeq <- 0  # Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
MassAtMaturityRelativeSeq <- 0.5  # Relative mass in relation to maximum Size
HeightBreadthSeq <- 1  # Relative height
MeanHeightSeq <- 0.5
DispersalKernelSeq <- 0  # The higher this values, the more local is the dispersal
DispersalKernelAsymmetrySeq <- 0.5  # The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)

# ============================================================================
# Define traits if neutral species pool(SpeciesPoolType=2) is choose
# The traits are the same for all species and have to be specified below. If correlations between traits are choosen,
# the age at maturity and the recruitment are based on the correlations defined above instead of the one defined below
MaxMassNeutral <- 100  # maximum mass of species/functional types (g)
AgeAtMaturityNeutral <- 3  # age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Neutral <- 15  # If no correlation is choose, this value is used as recruitment (it is not multiplied by RecruitmentInvestmentRel!)
RecruitmentInvestmentRelMeanNeutral <- 0.1
RecruitmentIncNeutral <- 0  # Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
MassAtMaturityRelativeNeutral <- 0.5  # Relative mass in relation to maximum Size
HeightBreadthNeutral <- 0.5  # Relative height
MeanHeightNeutral <- 0.5  # Mean height
DispersalKernelNeutral <- 1.5  # The higher this values, the more local is the dispersal
DispersalKernelAsymmetryNeutral <- 0.5  # The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)

