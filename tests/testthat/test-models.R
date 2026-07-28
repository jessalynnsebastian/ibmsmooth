test_that("only the two intended Stan implementations remain", {
  stan_dir <- system.file("stan", package = "ibmsmooth")
  if (!nzchar(stan_dir)) {
    stan_dir <- testthat::test_path("..", "..", "inst", "stan")
  }
  files <- sort(basename(list.files(stan_dir, pattern = "\\.stan$",
                                    full.names = TRUE)))
  expect_identical(files, c("ibm.stan", "ibm_adaptive.stan"))
})

test_that("ordinary IBM accepts an arbitrary Stan smoothing prior", {
  prior <- "tau ~ student_t(4, 0, 0.7);"
  code <- ibm(get_code = TRUE, smoothing_prior = prior)
  expect_match(code, prior, fixed = TRUE)
  expect_false(grepl("SMOOTHING_PRIOR", code, fixed = TRUE))
  expect_match(code, "real<lower=0> tau", fixed = TRUE)
})

test_that("adaptive Stan model matches the chapter hierarchy", {
  code <- ibm(get_code = TRUE, adaptive = TRUE)
  expect_match(code, "real tau_i = gamma * xi[i - 1]", fixed = TRUE)
  expect_match(code, "lambda_interval = square(gamma * xi)", fixed = TRUE)
  expect_match(
    code,
    "gamma = global_scale * tan(0.5 * pi() * gamma_unif)",
    fixed = TRUE
  )
  expect_match(code, "xi = tan(0.5 * pi() * xi_unif)", fixed = TRUE)
  expect_match(code, "gamma_unif ~ uniform(0, 1)", fixed = TRUE)
  expect_match(code, "xi_unif ~ uniform(0, 1)", fixed = TRUE)
  expect_false(grepl("~ cauchy", code, fixed = TRUE))
  expect_match(code, "h * fprime[i - 1]", fixed = TRUE)
  expect_match(code, "tau_i * sqrt_deltat[i - 1]", fixed = TRUE)
  expect_match(code, "generated quantities", fixed = TRUE)
  expect_match(code, "y_obs ~ normal(f[obs_time_idx], sigma)", fixed = TRUE)
  expect_false(grepl("operational", code, ignore.case = TRUE))
})

test_that("analytic Cholesky factor gives lambda Q(delta)", {
  for (delta in c(0.01, 0.4, 3.7)) {
    for (lambda in c(0.2, 2)) {
      L <- ibm_transition_cholesky(delta, lambda)
      expect_equal(
        L %*% t(L),
        ibm_transition_covariance(delta, lambda),
        tolerance = 1e-14
      )
    }
  }
})

test_that("adaptive simulation uses interval-specific exact transitions", {
  t <- c(0, 0.2, 1.5, 2)
  lambda <- c(0.1, 2, 0.3)
  sim <- simulate_ibm(t, lambda, initial = c(1, -1), seed = 4)
  expect_equal(sim$lambda_interval, lambda)
  expect_equal(unname(sim$state[1, ]), c(1, -1))
  expect_equal(sim$fprime, sim$state[, 1])
  expect_equal(sim$f, sim$state[, 2])
})

test_that("fit input validation is performed before sampling", {
  expect_error(ibm(1:3, c(1, 1, 1)), "positive finite standard deviation")
  expect_error(ibm(c(1, 1), c(0, 1)), "At least two unique")
  expect_error(ibm(get_code = TRUE, adaptive = NA), "adaptive")
  expect_error(
    ibm(get_code = TRUE, smoothing_prior = ""),
    "one non-empty character"
  )
})
