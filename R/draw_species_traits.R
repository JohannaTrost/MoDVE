#' Draw species traits from generating parameters
#'
#' @param species_params a list of trait-generating parameter, which must
#' contain the following elements:
#' \itemize{
#'  \item MaxMassLogScaleRandom TRUE/FALSE
#'  \item MaxMassRandom
#'  \item MassAtMaturityRelativeRandom vector of 2 elements
#'  \item InterceptAgeMaturity
#'  \item ScalingAgeMaturity
#'  \item AgeAtMaturityDevCorr
#'  \item CorrelationMassRecruitment TRUE/FALSE
#'  \item RecruitmentInvestmentRelMeanCorr
#'  \item RecruitmentInvestmentRelDevCorr
#'  \item RecruitmentNormalizeAtSize1Corr
#'  \item RecruitmentInvestmentRelMeanRandom vector of 2 elements
#'  \item RecruitmentNormalizeAtSize1Random vector of 2 elements
#'  \item DispersalKernelRandom
#'  \item DispersalKernelAsymmetryRandom
#'  \item HeightBreadthRandom
#'  \item HeightBreadthTheoretical
#'  \item Imax
#'  \item LAI
#'  \item kL
#' }
#'
#' @returns a named list of numeric containing the following traits:
#' \itemize{
#'  \item MaximumMass
#'  \item MassAtMaturity
#'  \item GrowthRate
#'  \item DispersalKernel
#'  \item DispersalKernelAsymmetry
#'  \item RecruitmentInvestmentRel
#'  \item RecruitmentInc
#'  \item MinLight
#'  \item MaxLight
#'  \item OptimumLight
#'  \item LightBreadth
#'  \item LightResponseA
#'  \item LightResponseB
#'  \item LightResponseC
#'  \item MinHeightRel
#'  \item MaxHeightRel
#'  \item MeanHeightRel,
#'  \item HeightBreadth
#'  \item MaxRecruitsAtMassAtMaturity
#'  \item AgeAtmaturity
#' }
#' @export
#'
draw_species_traits <- function(species_params) {

  # Unpack parameters
  check_species_params(species_params)
  list2env(config, envir = environment())

  # Draw max size
  if (MaxMassLogScaleRandom) {
    MaxMassLog <- runif(1, min = log10(MaxMassRandom[1]),
                        max = log10(MaxMassRandom[2]))
    MaxMass <- 10^MaxMassLog
  } else {
    MaxMass <- runif(1, min = MaxMassRandom[1], max = MaxMassRandom[2])
  }

  # Mass at maturity is a function of the maximum size
  MassAtMaturity <- MaxMass * runif(1, min = MassAtMaturityRelativeRandom[1],
                                    max = MassAtMaturityRelativeRandom[2])

  AgeAtMaturity <- InterceptAgeMaturity * (MaxMass^ScalingAgeMaturity) *
    runif(1, min = 1 - AgeAtMaturityDevCorr, max = 1 + AgeAtMaturityDevCorr)

  # Growth rate of the Bertalanffy growth curve
  K <- -(log(1) + log(1 - (MassAtMaturity / MaxMass))) / AgeAtMaturity

  # Recruitment
  if (CorrelationMassRecruitment) {
    RecruitmentInvestmentRel <- runif(1,
      RecruitmentInvestmentRelMeanCorr * (1 - RecruitmentInvestmentRelDevCorr),
      RecruitmentInvestmentRelMeanCorr * (1 + RecruitmentInvestmentRelDevCorr)
    )
    RecruitmentNormalizeAtSize1 <- RecruitmentNormalizeAtSize1Corr  # Factor converting the reproductive biomass to potential recruits
    SlopeRecruitment <- 0  # Slope of the correlation between mass and recruitment
    InterceptRecruitment <- RecruitmentNormalizeAtSize1
    RecruitmentInc <- 0
  } else {
    RecruitmentInvestmentRel <- runif(1,
      RecruitmentInvestmentRelMeanRandom[1],
      RecruitmentInvestmentRelMeanRandom[2]
    )
    RecruitmentNormalizeAtSize1 <- runif(1,
      RecruitmentNormalizeAtSize1Random[1],
      RecruitmentNormalizeAtSize1Random[2]
    )
    SlopeRecruitment <- 0  # No slope if no correlation is choosen
    InterceptRecruitment <- RecruitmentNormalizeAtSize1 - SlopeRecruitment
    RecruitmentInc <- runif(1, RecruitmentIncRandom[1], RecruitmentIncRandom[2])
      # Not meaningful if no correlation
  }
  MaxRecruits <- InterceptRecruitment * RecruitmentInvestmentRel
  MaxRecruitsMaturity <- (InterceptRecruitment + SlopeRecruitment *
                            MassAtMaturity) * RecruitmentInvestmentRel
  # Dispersal
  DispersalKernel <- runif(1, DispersalKernelRandom[1], DispersalKernelRandom[2])
  DispersalKernelAsymmetry <- runif(1, DispersalKernelAsymmetryRandom[1], DispersalKernelAsymmetryRandom[2])

  # Height niche
  MeanHeight <- runif(1, min = 0, max = 1)  # relative height in relation to canopy height
  HeightBreadthTheoretical <- runif(1, HeightBreadthRandom[1], HeightBreadthRandom[2])
  MinHeight <- max(c(0, MeanHeight - (HeightBreadthTheoretical / 2)))
  MaxHeight <- min(c(1, MeanHeight + (HeightBreadthTheoretical / 2)))
  HeightBreadth <- MaxHeight - MinHeight

  # Light niche
  MinLight <- Imax * exp(-kL * LAI * (1 - MinHeight))
  MaxLight <- Imax * exp(-kL * LAI * (1 - MaxHeight))
  OptimumLight <- (MaxLight + MinLight) / 2
  LightBreadth <- MaxLight - MinLight
  light_resp_params <- get_light_resp_params(MinLight, MaxLight, OptimumLight)

  # Output
  sp_traits <- list(
    MaxMass, MassAtMaturity, K, DispersalKernel, DispersalKernelAsymmetry,
    RecruitmentInvestmentRel, RecruitmentInc, MinLight, MaxLight, OptimumLight,
    LightBreadth, light_resp_params[1], light_resp_params[2], light_resp_params[3],
    MinHeight, MaxHeight, MeanHeight, HeightBreadth,
    MaxRecruits, MaxRecruitsMaturity, AgeAtMaturity
  )
  names(sp_traits) <- species_trait_names()
  return(sp_traits)
}
