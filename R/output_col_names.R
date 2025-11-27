species_trait_names <- function() {
  return(c("MaximumMass", "MassAtMaturity", "GrowthRate",
           "DispersalKernel", "DispersalKernelAsymmetry",
           "RecruitmentInvestmentRel", "RecruitmentInc",
           "MinLight", "MaxLight", "OptimumLight", "LightBreadth",
           "LightResponseA", "LightResponseB", "LightResponseC",
           "MinHeightRel",  "MaxHeightRel", "MeanHeightRel",
           "HeightBreadth", "MaxRecruitsAtMaxMass",
           "MaxRecruitsAtMassAtMaturity", "AgeAtMaturity"))
}

inds_output_names <- function() {
  return(c(
    "SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z",
    "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel"
    ))
}

species_output_names <- function() {
  return(c(
    "TimeStep", "SpeciesID", "NumberIndividualsBeginning",
    "NumberIndividualsEnd", "NumberMatureIndividuals", "NumberRecruits",
    "NumberRecruitsPotential", "NumberMortalityBranchFall",
    "NumberMortalityLight", "NumberMortalityCompetition",
    "NumberMortalityNatural", "PopulationGrowthRate",
    "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass",
    "AverageAge", "MinLight", "MaxLight",
    "MeanLight", "MinHeight", "MaxHeight", "MeanHeight"
  ))
}

comm_output_names <- function() {
  return(c(
    "timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
    "NumberIndividualsBeginning", "NumberIndividualsEnd", "nb_recruits_matrix",
    "MortalityBranchFall", "MortalityLight", "MortalityCompetition",
    "MortalityNatural", "BranchSurfaceIndex", "EpiphyteFilling"
  ))
}
