#' Laplace Transform of IV
#'
#' @description
#' Laplace transform of Integrated Variance (IV) with condition on two side
#' levels \eqn{v_u} and \eqn{v_t}. It is different from Moment-Generating
#' function:
#' \eqn{
#'   \mathcal{L}(a) = M(-a).
#' }
#' Closely related with the Characteristic Function through
#' \eqn{\Phi(a) = \mathcal{L}(-ia)}, where \eqn{\Phi(a)} is the Characteristic
#' Function which is implemented as the function [ajd.sim.bk::CF()].
#'
#' @param a Laplace transform input, can be either a real number or a complex
#' number.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#'
#' @return Laplace transform result, either a scalar real number or a scalar
#' complex number, corresponding to the input whether a real number or a
#' complex number.
#' @export
#'
#' @seealso Closely related with the Characteristic Function [ajd.sim.bk::CF()].
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' a = 1
#' Laplace(a, v0, v1, tau, k, theta, sigma)
Laplace <- function(a, v0, v1, tau, k, theta, sigma) {
  # if (is.complex(a)) {
  #   cat(sprintf("\na = (%.7f, %.7fi)\n", Re(a), Im(a)))
  #   x = k^2 + 2*sigma^2*a
  #   cat(sprintf("k^2 + 2*sigma^2*a = (%.7f, %.7fi)\n", Re(x), Im(x)))
  #   cat(sprintf("sqrt(k^2 + 2*sigma^2*a) = (%.7f, %.7fi)\n", Re(sqrt(x)), Im(sqrt(x))))
  # } else {
  #   cat(sprintf("\na = %.7f, k = %.7f, sigma = %.7f\n", a, k, sigma))
  #   x = k^2 + 2*sigma^2*a
  #   cat(sprintf("k^2 + 2*sigma^2*a = %.7f\n", x))
  # }
  a = ifelse(is.complex(a), a, as.complex(a))
  #
  g_a = sqrt(k^2 + 2*sigma^2*a) # gamma(a)
  dk = exp(-k*tau); da = exp(-g_a*tau)
  #
  f = (g_a/k) * exp(-0.5*(g_a-k)*tau) * ((1-dk)/(1-da))
  #
  t1 = (v0+v1)/sigma^2
  t2 = k * ((1+dk)/(1-dk)) - g_a * ((1+da)/(1-da))
  f = f * exp(t1 * t2)
  #
  nu = 2*theta*k/sigma^2 - 1
  x1 = sqrt(v0*v1) * (4*g_a/sigma^2) * exp(-0.5*g_a*tau)/(1-da)
  x2  = sqrt(v0*v1) * (4*k  /sigma^2) * exp(-0.5*k*tau) / (1-dk)
  #
  if (is.complex(x1)) {
    num_scaled = abs(Re(x1)) > 400; cond1 = abs(Re(x1)) < 1000
  } else {
    num_scaled = abs(x1) > 400;     cond1 = abs(x1) < 1000
  }
  mybesselI = ifelse(cond1, Bessel::BesselI, Bessel::besselIasym)
  num = mybesselI(x1, nu, expon.scaled = num_scaled)
  #
  den_scaled = abs(x2) > 400
  mybesselI = ifelse(abs(x2) < 1000, Bessel::BesselI, Bessel::besselIasym)
  den = mybesselI(x2, nu, expon.scaled = den_scaled)
  f = f * (num/den)
  #
  if (is.complex(x1)) {
    f = f * exp(num_scaled * abs(Re(x1)) - den_scaled * abs(x2))
  } else {
    f = f * exp(num_scaled * abs(x1) - den_scaled * abs(x2))
  }
  # f = f * (BesselI(as.complex(x1), nu) / besselI(x2, nu))
  return(f)
}
