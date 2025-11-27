#' Title
#'
#' @param E
#' @param Microhabitat
#' @param MortalityMethod
#' @param MortRateRandom
#' @param MortRateMass
#' @param MortRateMassScaling
#'
#' @returns
#' @export
#'
#' @examples
resolve_mortality <- function(E, Microhabitat, MortalityMethod, MortRateRandom,
                  MortRateMass, MortRateMassScaling) {

  for (i in seq_len(nrow(E))) {
    if (E$Status[i] == 1) {

      vox <- Microhabitat[E$X[i], E$Y[i], E$Z[i],]

      # The following comparison would fail without the is.nan check,
      # because Microhabitat contains NaNs in some entries and
      # in R a comparison with a NaN returns NA, not a boolean.
      # Note: We call runif repeatedly intentionally. See Issue #16 on Github
      # TODO: discuss priority among the different sources of mortality
      # e.g. the rate of random mortality won't match the parameter
      # because some fraction has already died from branch fall/light
      # Once priority is clarified drawing mortality from sa loss and
      # random/mass mortality could be vectorised

      # Branch fall mortality
      if (!is.nan(vox[2]) && runif(1, min = 0, max = 1) < vox[2]) {
        E$Status[i] <- 3
      } else if (vox[3] < E$MinLight[i] | vox[3] > E$MaxLight[i]) {
        # Unsuitable light conditions
        E$Status[i] <- 4
      } else if (MortalityMethod == 0 && runif(1, min = 0, max = 1) < MortRateRandom) {  # Natural mortality rate
        # Baseline random mortality
        E$Status[i] <- 5
      } else if (MortalityMethod == 1 && runif(1, min = 0, max = 1) < (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
        # Mass-dependent mortality
        E$Status[i] <- 5
      }
    }
  }
  return(E)
}
