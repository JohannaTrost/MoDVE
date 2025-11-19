#' Derive coefficients of the growth-light response from light parameters
#'
#' The growth response is a parabolic function of light intensity:
#' \deqn{y = a * light^2 + b*light + c}
#'
#' This function derives the coefficients
#'
#' @param MinLight numeric, the minimum light in which the species can survive
#' @param MaxLight numeric, the maximum light at which the species can suvive
#' @param OptimumLight numeric, light intensity at which growth is maximized
#'
#' @returns a vector of three numerics, parameters a, b, c of the parabolic function
#' @export
#'
get_light_resp_params <- function(MinLight, MaxLight, OptimumLight) {

  # Parameters of parabolic response curve y = ax^2 + bx + c such that:
  # f(MinLight) = 0
  # f(MaxLight) = 0
  # f(OptimumLight) = 1

  x1 <- MinLight
  y1 <- 0

  x2 <- MaxLight
  y2 <- 0

  x3 <- OptimumLight
  y3 <- 1

  # Derive coefficients from the three known points
  # see Petter et al. 2021 Appendix A2
  a <- (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) /
    ((x1 - x2) * (x1 - x3) * (x3 - x2)) # eq 24

  b <- (x1^2 * (y2 - y3) + x2^2 * (y3 - y1) + x3^2 * (y1 - y2)) /
    ((x1 - x2) * (x1 - x3) * (x2 - x3)) # eq 25

  c <- (x1^2 * (x2 * y3 - x3 * y2) +
          x1 * (x3^2 * y2 - x2^2 * y3) +
          x2 * x3 * y1 * (x2 - x3)
  ) / ((x1 - x2) * (x1 - x3) * (x2 - x3)) # eq 26

  return(c(a, b, c))
}
