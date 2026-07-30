test_that("non-centered and centered Stan implementations are present", {
  stan_dir <- system.file("stan", package = "ibmsmooth")
  if (!nzchar(stan_dir)) {
    stan_dir <- testthat::test_path("..", "..", "inst", "stan")
  }
  files <- sort(basename(list.files(stan_dir, pattern = "\\.stan$",
                                    full.names = TRUE)))
  expect_identical(
    files,
    c("ibm.stan", "ibm_adaptive.stan",
      "ibm_adaptive_centered.stan",
      "ibm_adaptive_centered_regularized.stan",
      "ibm_adaptive_regularized.stan", "ibm_centered.stan")
  )
})

test_that("horseshoe is default and regularized horseshoe is selectable", {
  default_code <- ibm(get_code = TRUE, adaptive = TRUE)
  regularized_code <- ibm(
    get_code = TRUE, adaptive = TRUE,
    adaptive_prior = "regularized_horseshoe"
  )
  centered_code <- ibm(
    get_code = TRUE, adaptive = TRUE,
    adaptive_prior = "regularized_horseshoe",
    parameterization = "centered"
  )
  expect_false(grepl("slab_aux", default_code, fixed = TRUE))
  for (code in list(regularized_code, centered_code)) {
    expect_match(code, "slab_aux ~ inv_gamma", fixed = TRUE)
    expect_match(code, "xi_regularized", fixed = TRUE)
    expect_match(
      code,
      "square(slab) + square(gamma) * square(xi)",
      fixed = TRUE
    )
    expect_match(code, "square(gamma * xi_regularized)", fixed = TRUE)
  }
})

test_that("regularized retry is explicit and non-applicable fits are ignored", {
  expect_identical(
    eval(formals(ibm)$regularized_retry),
    c("ask", "never")
  )
  expect_error(
    ibm(get_code = TRUE, regularized_retry = "always"),
    "arg"
  )
  mock <- structure(
    list(stanfit = list(), data = list(), adaptive = FALSE),
    class = c("ibmfit", "list")
  )
  check <- check_regularized_horseshoe(mock)
  expect_false(check$recommended)
})

test_that("centered programs use the same exact state transition", {
  for (adaptive in c(FALSE, TRUE)) {
    code <- ibm(
      get_code = TRUE, adaptive = adaptive, parameterization = "centered"
    )
    expect_match(code, "vector[T] fprime", fixed = TRUE)
    expect_match(code, "vector[T] f", fixed = TRUE)
    expect_match(code, "multi_normal_cholesky", fixed = TRUE)
    expect_match(code, "h * sqrt(h) / 2", fixed = TRUE)
    expect_match(code, "inv_sqrt(12.0)", fixed = TRUE)
    expect_match(code, "vector[T] fprime", fixed = TRUE)
    expect_false(grepl("array[T - 1] vector[2] z_transition",
                       code, fixed = TRUE))
  }
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
  expect_match(code, "global_prior == 1", fixed = TRUE)
  expect_match(code, "global_scale * tan", fixed = TRUE)
  expect_match(code, "global_scale * inv_Phi", fixed = TRUE)
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

test_that("reference-process global priors are calibrated correctly", {
  grid <- c(0, 0.1, 0.4, 1)
  reference_sd <- ibmsmooth:::.ibm_reference_sd(grid)
  expect_equal(
    reference_sd,
    exp(mean(log(sqrt(grid[-1]^3 / 3))))
  )

  upper <- 1.2
  alpha <- 0.1
  normal_scale <- ibmsmooth:::.ibm_global_scale(
    reference_sd, upper, alpha, "half_normal"
  )
  cauchy_scale <- ibmsmooth:::.ibm_global_scale(
    reference_sd, upper, alpha, "half_cauchy"
  )
  expect_equal(
    2 * stats::pnorm(-upper / (reference_sd * normal_scale)),
    alpha
  )
  expect_equal(
    1 - 2 / pi * atan(upper / (reference_sd * cauchy_scale)),
    alpha
  )
})

test_that("all Stan models can omit the likelihood for prior sampling", {
  for (adaptive in c(FALSE, TRUE)) {
    for (parameterization in c("noncentered", "centered")) {
      code <- ibm(
        get_code = TRUE, adaptive = adaptive,
        parameterization = parameterization
      )
      expect_match(code, "int<lower=0, upper=1> prior_only", fixed = TRUE)
      expect_match(code, "if (!prior_only)", fixed = TRUE)
    }
  }
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

test_that("IBM bridge covariance has the expected endpoints", {
  h <- 1.7
  for (u in c(0, h)) {
    cross <- ibmsmooth:::.ibm_Q(u) %*% t(ibmsmooth:::.ibm_A(h - u))
    gain <- cross %*% solve(ibmsmooth:::.ibm_Q(h))
    conditional <- ibmsmooth:::.ibm_Q(u) - gain %*% t(cross)
    expect_equal(conditional, matrix(0, 2, 2), tolerance = 1e-12)
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
    ibm(get_code = TRUE, parameterization = "hybrid"),
    "arg"
  )
  expect_error(
    ibm(get_code = TRUE, adaptive = TRUE, adaptive_prior = "other"),
    "arg"
  )
  expect_error(
    ibm(get_code = TRUE, adaptive = TRUE, global_prior = "other"),
    "arg"
  )
  expect_error(
    ibm(get_code = TRUE, smoothing_prior = ""),
    "one non-empty character"
  )
})
