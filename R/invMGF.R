#' Compute the Moments of IV
#'
#' @description
#' Compute the moments of IV by numerical inversion of an adaptively modified
#' moment generating function introduced by Choudhury and Lucantoni (1996).
#'
#' @param m largest order of the moments to be computed.
#' @param gamma accuracy parameter, error bound \eqn{10^{-\gamma}}.
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @return vector of moments, \eqn{\mu_1,\cdots,\mu_m}.
#' @export
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' m = 4; gamma = 11
#' mu = invMGF(m, gamma, v0, v1, tau, k, theta, sigma)
invMGF <- function(m, gamma, v0, v1, tau, k, theta, sigma) {
  mu = rep(0, m)
  #
  l = 1; n = 1; a_n = 1; r_n = 10^(-gamma/(2*n*l))
  mu[1] = factorial(n)/(2*n*l * r_n^n * a_n^n)
  cf1 = Laplace(-a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf2 = Laplace( a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf3 = 0
  mu[1] = mu[1] * (cf1 + (-1)^n*cf2 + 2*cf3)
  # print(sprintf("mu[1]=%.10f", mu[1]))
  #
  n = 2; a_n = 1/mu[1]; r_n = 10^(-gamma/(2*n*l))
  mu[2] = factorial(n)/(2*n*l * r_n^n * a_n^n)
  cf1 = Laplace(-a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf2 = Laplace( a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf3 = 0
  for (j in 1:(n*l -1)) {
    z = Laplace(-a_n*r_n * exp(1i*pi*j/(n*l)), v0, v1, tau, k, theta, sigma)
    z = z * exp(-1i*pi*j/l)
    cf3 = cf3 + Re(z)
  }
  mu[2] = mu[2] * (cf1 + (-1)^n*cf2 + 2*cf3)
  #
  # update mu1 and mu2
  #
  # l = max(1,2);
  l = 1; a_n = 2*mu[1]/mu[2]
  #
  n = 1; r_n = 10^(-gamma/(2*n*l))
  mu[1] = factorial(n)/(2*n*l * r_n^n * a_n^n)
  cf1 = Laplace(-a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf2 = Laplace( a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf3 = 0
  mu[1] = mu[1] * (cf1 + (-1)^n*cf2 + 2*cf3)
  # print(sprintf("mu[1]=%.10f", mu[1]))
  #
  n = 2; r_n = 10^(-gamma/(2*n*l))
  mu[2] = factorial(n)/(2*n*l * r_n^n * a_n^n)
  cf1 = Laplace(-a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf2 = Laplace( a_n*r_n, v0, v1, tau, k, theta, sigma)
  cf3 = 0
  for (j in 1:(n*l - 1)) {
    z = Laplace(-a_n*r_n * exp(1i*pi*j/(n*l)), v0, v1, tau, k, theta, sigma)
    z = z * exp(-1i*pi*j/l)
    cf3 = cf3 + Re(z)
  }
  mu[2] = mu[2] * (cf1 + (-1)^n*cf2 + 2*cf3)
  #
  if (m < 3) { stop(sprintf("m must >=3, however supplied with m=%i", m)) }
  for (n in 3:m) {
    # l = max(1, 2);
    l = 1; a_n = (n-1)*mu[n-2]/mu[n-1]; r_n = 10^(-gamma/(2*n*l))
    mu[n] = factorial(n)/(2*n*l * r_n^n * a_n^n)
    cf1 = Laplace(-a_n*r_n, v0, v1, tau, k, theta, sigma)
    cf2 = Laplace( a_n*r_n, v0, v1, tau, k, theta, sigma)
    cf3 = 0
    for (j in 1:(n*l - 1)) {
      z = Laplace(-a_n*r_n * exp(1i*pi*j/(n*l)), v0, v1, tau, k, theta, sigma)
      z = z * exp(-1i*pi*j/l)
      cf3 = cf3 + Re(z)
    }
    mu[n] = mu[n] * (cf1 + (-1)^n*cf2 + 2*cf3)
  }
  #
  mu_new = rep(0, length(mu))
  for (i in 1:length(mu)) {
    if (Im(mu[i]) > 1e-5) stop(sprintf("Im(mu[%d]) > 1e-5", Im(mu[i])))
    mu_new[i] = Re(mu[i])
  }
  return(mu_new)
}
