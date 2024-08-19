#' Simulating IV from Conditional Distribution
#'
#' @description
#' Random variable generation for the conditional distribution of the
#' Integrated Variance (IV), through a Pearson Distribution approximation
#' which match the first four moments of the unknown true conditional
#' distribution of IV.
#'
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @details
#' The first four moments are calculated through [invMGF()].
#' **Note that:**
#' - when the time difference `tau` is very short, like < 0.0001, it is very
#'   likely to introduce computational errors, therefore, we let
#'   `iv = tau * (v0+v1)/2` directly.
#'
#' - when the computed variance is very small, the computed skewness and
#'   kurtosis are highly un-reliable, therefore, we let `iv` equal the mean
#'   directly.
#'
#' @return scalar
#' @export
#'
#' @seealso See [ajd.sim.bk::riv()] for the Broadie-Kaya algorithm which
#' is through Fourier inversion of its Characteristic Function.
#'
#' @references
#' 1. Kyriakou, I., Brignone, R., & Fusai, G. (2024). Unified moment-based modeling
#'   of integrated stochastic processes. *Operations Research*, 72(4), 1630-1653.
#'
#' 2. Broadie, M., & Kaya, Ö. (2006). Exact simulation of stochastic volatility
#'   and other affine jump diffusion processes. *Operations Research*,
#'   54(2), 217-231.
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' iv = riv(v0, v1, tau, k, theta, sigma)
riv <- function(v0, v1, tau, k, theta, sigma) {
  if (tau < 0.0001) { return(tau * (v0+v1)/2) }
  #
  mu = invMGF(4, 11, v0, v1, tau, k, theta, sigma)
  stdmoment = stdmom(mu)
  #
  mean = stdmoment[1]; var = stdmoment[2]
  skew = stdmoment[3]; kurt = stdmoment[4]
  #
  if (var < 1e-5) { return(mean) }
  #
  # check kurtosis bound
  if (kurt < skew^2 + 1) {
    temp = "\nriv(%.7f, %.7f, %.7f, %.4f, %.4f, %.4f)\n"
    cat(sprintf(temp, v0, v1, tau, k, theta, sigma))
    #
    temp = "mean=%.7f, var=%.20f, skew=%.7f, kurt=%.7f\n"
    cat(sprintf(temp, mean, var, skew, kurt))
    #
    temp = "increase kurtosis from (=%.10f) to its bound (=%.10f)\n"
    cat(sprintf(temp, kurt, skew^2+1))
    # adjustment
    kurt = skew^2 + 1; stdmoment[4] = kurt
  }
  # sample IV from Pearson Distribution
  iv = PearsonDS::rpearson(1, moments = stdmoment)
  #
  return(iv)
}
