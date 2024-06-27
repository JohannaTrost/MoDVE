# Create species matrices

AgeMaturityMetabolic <- function(InterceptAgeMaturity, ScalingAgeMaturity, Mass) {
    return(InterceptAgeMaturity * (Mass^ScalingAgeMaturity))
}

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

# ============================================================================
# Create folder to save species trait matrices
if (SpeciesPoolType == 0) {
    if (CorrelationMassAgeOfMaturity == 1 || CorrelationMassRecruitment == 1) {
        FullNameSpeciesPool <- paste("SP_Random_", NameSpeciesPool, "_TraitCorrOn", sep="")
    } else {
        FullNameSpeciesPool <- paste("SP_Random_", NameSpeciesPool, "_TraitCorrOff", sep="")
    }
} else if (SpeciesPoolType == 1) {
    if (CorrelationMassAgeOfMaturity == 1 || CorrelationMassRecruitment == 1) {
        FullNameSpeciesPool <- paste("SP_Sequential_", NameSpeciesPool, "_TraitCorrOn", sep="")
    } else {
        FullNameSpeciesPool <- paste("SP_Sequential_", NameSpeciesPool, "_TraitCorrOff", sep="")
    }
} else if (SpeciesPoolType == 2) {
    if (CorrelationMassAgeOfMaturity == 1 || CorrelationMassRecruitment == 1) {
        FullNameSpeciesPool <- paste("SP_Neutral_", NameSpeciesPool, "_TraitCorrOn", sep="")
    } else {
        FullNameSpeciesPool <- paste("SP_Neutral_", NameSpeciesPool, "_TraitCorrOff", sep="")
    }
}

# where to save
SaveDirectory <- file.path(MainOutputDirectory, FullNameSpeciesPool)
dir.create(SaveDirectory, recursive=TRUE)

ColumnHeaders <- c("SpeciesID", "MaximumMass", "MassAtMaturity", "GrowthRate",
                   "DispersalKernel", "DispersalKernelAsymmetry", "RecruitmentInvestmentRel",
                   "RecruitmentInc", "MinLight", "MaxLight", "OptimumLight", "LightBreadth",
                   "LightResponseA", "LightResponseB", "LightResponseC", "MinHeightRel",
                   "MaxHeightRel", "MeanHeightRel", "HeightBreadth", "MaxRecruitsAtMaxMass",
                   "MaxRecruitsAtMassAtMaturity", "AgeAtMaturity")

# Main loop (for random generation of species pool)
for (Num in 1:numSpeciesPools) {

    # Trait matrix where the trait information of each species is saved
    SpeciesTraitMatrix <- array(rep(0, NumberOfSpecies * 16), dim=c(NumberOfSpecies, 16))

    for (NumSpecies in 1:NumberOfSpecies) {

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
        RecruitmentNormalizeAtSize1Corr <- 70  # Factor converting the reproductive biomass to potential recruits
        SlopeRecruitmentCorr <- 0 * RecruitmentNormalizeAtSize1Corr  # Slope of the correlation between mass and recruitment
        RecruitmentInvestmentRelMeanCorr <- 0.1  # Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
        RecruitmentInvestmentRelDevCorr <- 0.25  # The relative deviation from the mean recruitment
        RecruitmentIncMaxCorr <- 0

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

        DispersalKernel <- runif(1, min=DispersalKernelRandom[1], max=DispersalKernelRandom[2])
        DispersalKernelAsymmetry <- runif(1, min=DispersalKernelAsymmetryRandom[1], max=DispersalKernelAsymmetryRandom[2])

        # ============================================================================
        # Traits of ecologcial niche
        # 1. Randomly choose mean height and height breadth
        MeanHeight <- runif(1, min=0, max=1)  # realtive height in relation to canopy height
        HeightBreadthTheoretical <- runif(1, min=HeightBreadthRandom[1], max=HeightBreadthRandom[2])

        # Minimum and maximum height under which the species is
        # able to survive
        MinHeight <- max(c(0, MeanHeight - (HeightBreadthTheoretical / 2)))
        MaxHeight <- min(c(1, MeanHeight + (HeightBreadthTheoretical / 2)))
        HeightBreadth <- MaxHeight - MinHeight

        # 2. Convert heigth ranges to light ranges
        # For this, a standard forest with the following parameters is
        # assumed
        MinLight <- Imax * exp(-kL * LAI * (1 - MinHeight))
        MaxLight <- Imax * exp(-kL * LAI * (1 - MaxHeight))
        OptimumLight <- (MaxLight + MinLight) / 2
        LightBreadth <- MaxLight - MinLight

        # 3. Calculate parameters of parabolic response curve y=ax^2+bx+c
        # We assume that the function is a paraboloid which goes trough
        # three points (MinLight/0) (MaxLight/0) (OptimumLight/1)

        x1 <- MinLight
        x2 <- MaxLight
        x3 <- OptimumLight

        y1 <- 0
        y2 <- 0
        y3 <- 1

        a <- (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / ((x1-x2) * (x1-x3) * (x3-x2))
        b <- (x1^2*(y2-y3) + x2^2*(y3-y1) + x3^2*(y1-y2)) / ((x1-x2) * (x1-x3) * (x2-x3))
        c <- (x1^2*(x2*y3 - x3*y2) + x1*(x3^2*y2 - x2^2*y3) + x2*x3*y1*(x2-x3)) / ((x1-x2) * (x1-x3) * (x2-x3))
        # In the model, based on these parameters the light response for each species can be calculated:
        # Parabol=@(a,b,c,x) a*x^2+b*x+c;

}
