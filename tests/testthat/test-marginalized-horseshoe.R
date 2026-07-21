dense_state_posterior <- function(y, obs_time_idx, dt, sigma, tau,
                                  initial_sd = 5) {
  TT <- length(dt) + 1L
  G <- matrix(0, 2L * TT, 2L * TT)
  G[1:2, 1:2] <- initial_sd^2 *
    matrix(c(1, 0.5, 0.5, 1 / 3), 2L, 2L)
  for (i in 2:TT) {
    h <- dt[i - 1L]
    F <- matrix(c(1, h, 0, 1), 2L, 2L)
    S <- tau[i - 1L]^2 *
      matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2L, 2L)
    old <- seq_len(2L * (i - 1L))
    prev <- tail(old, 2L)
    now <- (2L * i - 1L):(2L * i)
    G[now, old] <- F %*% G[prev, old]
    G[old, now] <- t(G[now, old])
    G[now, now] <- F %*% G[prev, prev] %*% t(F) + S
  }
  H <- matrix(0, length(y), 2L * TT)
  H[cbind(seq_along(y), 2L * obs_time_idx)] <- 1
  V <- H %*% G %*% t(H) + diag(sigma^2, length(y))
  K <- G %*% t(H) %*% solve(V)
  list(
    mean = as.numeric(K %*% y),
    covariance = (G - K %*% H %*% G + t(G - K %*% H %*% G)) / 2
  )
}

test_that("Kalman likelihood equals dense Gaussian likelihood", {
  set.seed(10)
  t <- sort(runif(9))
  y <- rnorm(9)
  tau <- exp(rnorm(8, -2, 0.4))
  ans <- validate_marginalized_horseshoe(t, y, sigma = 0.3, tau = tau)
  expect_lt(abs(ans[["difference"]]), 1e-9)
})

test_that("simulation smoother matches dense moments with replicates and infer_at", {
  # Grid points 2 and 5 have no observations; point 3 has two replicates.
  dt <- c(0.6, 1.1, 0.4, 1.3)
  obs_time_idx <- c(1L, 3L, 3L, 4L)
  y <- c(-0.2, 0.7, 0.5, 1.1)
  sigma <- 0.25
  tau <- c(0.16, 0.09, 0.22, 0.12)
  exact <- dense_state_posterior(y, obs_time_idx, dt, sigma, tau)

  set.seed(14)
  draws <- replicate(
    6000,
    as.vector(ibmsmooth:::.ibmsmooth_ffbs(
      y, obs_time_idx, dt, sigma, tau, initial_sd = 5
    ))
  )
  empirical_mean <- rowMeans(draws)
  empirical_var <- apply(draws, 1L, stats::var)

  expect_equal(empirical_mean, exact$mean, tolerance = 0.035)
  expect_equal(empirical_var, diag(exact$covariance), tolerance = 0.025)
})

test_that("marginalized data preparation sorts observations and retains grid points", {
  prep <- ibmsmooth:::.ibmsmooth_marginalized_data(
    t = c(3, 1, 3, 5), y = c(2, 0, 1, 4), infer_at = c(2, 4),
    log_sigma = list(mu = -1, sd = 1), zeta = 0.1
  )
  expect_equal(prep$data$time_grid_raw, 1:5)
  expect_true(all(diff(prep$stan_data$obs_time_idx) >= 0))
  expect_equal(prep$stan_data$T, 5L)
  expect_equal(length(prep$stan_data$deltat), 4L)
})

test_that("ibm exposes the marginalized horseshoe as an explicit opt-in", {
  body_text <- paste(deparse(body(ibm)), collapse = "\n")
  expect_match(body_text, "ibm_horseshoe_marginalized", fixed = TRUE)
  expect_identical(eval(formals(ibm)$stan_horseshoe_engine),
                   c("joint", "marginalized"))
})
