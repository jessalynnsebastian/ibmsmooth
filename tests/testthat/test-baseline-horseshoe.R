test_that("baseline horseshoe is a separate selectable Stan model", {
  code <- ibm_smooth(adaptive = "baseline_horseshoe", get_code = TRUE,
                     centered = FALSE)
  expect_match(code, "square(tau0) + square(psi[i])", fixed = TRUE)
  expect_match(code, "tau[i - 1] * Lt[i - 1]", fixed = TRUE)

  original <- ibm_smooth(adaptive = "horseshoe", get_code = TRUE,
                         centered = FALSE)
  expect_false(grepl("tau0", original, fixed = TRUE))
})

test_that("baseline horseshoe has the requested variance limits", {
  dt <- c(0.2, 1.7, 3.1)
  tau0 <- 0.4
  psi <- c(0, 0.3, 1.2)
  ibm_covariance <- function(h) {
    matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2, 2)
  }
  covariance <- Map(function(h, p) (tau0^2 + p^2) * ibm_covariance(h),
                    dt, psi)

  expect_equal(covariance[[1]], tau0^2 * ibm_covariance(dt[1]))
  expect_equal(covariance[[2]] / (tau0^2 + psi[2]^2), ibm_covariance(dt[2]))
  expect_true(all(diag(covariance[[3]]) > diag(tau0^2 * ibm_covariance(dt[3]))))
})

test_that("centered baseline horseshoe preserves the transition mean", {
  code <- ibm_smooth(adaptive = "baseline_horseshoe", get_code = TRUE,
                     centered = TRUE)
  expect_match(code, "mu[1] = fprime[i - 1]", fixed = TRUE)
  expect_match(code,
               "mu[2] = f[i - 1] + deltat[i - 1] * fprime[i - 1]",
               fixed = TRUE)
  expect_match(code, "vector<lower=0, upper=1>[T-1] z_lambda", fixed = TRUE)
})
