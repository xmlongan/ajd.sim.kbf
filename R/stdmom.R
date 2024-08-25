#' Mean, Variance, Skewness and Kurtosis
#'
#' @description
#' Get the mean, variance, skewness and kurtosis from the first four
#' non-central moments. It should be noted that
#' - \eqn{variance > 0,} and
#' - \eqn{kurtosis \le skewness^2 + 1.}
#' When either of these two scenarios happen, a error is issued or other
#' modifications are conducted. When \eqn{std < 1e^{-5}}, a warning is issued,
#' warning that the skewness and kurtosis are highly un-reliable.
#'
#' @param mu vector of the first four non-central moments.
#'
#' @return vector of mean, variance, skewness and kurtosis.
#' @export
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' m = 4; gamma = 11
#' mu = invMGF(m, gamma, v0, v1, tau, k, theta, sigma)
#' stdmoment = stdmom(mu)
stdmom <- function(mu) {
  m1 = mu[1]; m2 = mu[2]; m3 = mu[3]; m4 = mu[4]
  variance = m2 - m1^2
  std = sqrt(variance)
  skewness = (m3 - 3*m1*m2 + 2*m1^3)/(std^3)
  kurtosis = (m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4)/(std^4)
  # check kurtosis bound and variance > 0
  if (variance <= 0) {
    stop(sprintf("variance(=%.7f) <= 0\n", variance))
  }
  if (std < 1e-5) {
    warning("std is < 1e-5, consequently,
            the skewness and kurtosis are highly un-reliable")
  }
  if (kurtosis < skewness^2 + 1) {
    cat(sprintf("std = %.7f\n", std))
    temp = "kurtosis (=%.7f) < skewness^2 (=%.7f) + 1\n"
    cat(sprintf(temp, kurtosis, skewness^2))
    temp = "I manually increase kurtosis from (=%.10f) to its bound (=%.10f)\n"
    cat(sprintf(temp, kurtosis, skewness^2+1))
    kurtosis = skewness^2 + 1
  }
  return(c(m1, variance, skewness, kurtosis))
}
