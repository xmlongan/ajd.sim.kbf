test_that("ryield_SVCJ works for most scenarios", {
  # expect_equal(2 * 2, 4)
  v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14; rho = -0.82
  r = 0.0319; tau = 1; lambda = 0.47; mu_bar = -0.1; sigma_s = 0.0001
  mu_v = 0.05; rho_J = -0.38
  N = 1000 # number of samples
  expect_no_error(ryield_SVCJ(N, v0, tau, r, k, theta, sigma, rho,
                              lambda, mu_bar, sigma_s, mu_v, rho_J))
})
