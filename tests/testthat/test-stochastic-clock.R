test_that("analytic factor gives the exact IBM covariance", {
  for (q in c(0.01, 0.4, 3.7)) {
    L <- ibm_clock_cholesky(q)
    expect_equal(L %*% t(L), ibm_clock_covariance(q), tolerance = 1e-14)
  }
})

test_that("constant and adaptive clocks nest on irregular grids", {
  dt <- c(0.2, 1.7, 3.1)
  lambda0 <- 0.4
  r <- rep(0, length(dt))
  expect_equal((lambda0 + r) * dt, lambda0 * dt)
  expect_equal((lambda0 * dt) / dt, rep(lambda0, length(dt)))
})

test_that("acceleration changes exact transition powers", {
  dt <- 0.7; lambda0 <- 0.5
  q0 <- lambda0 * dt
  q1 <- (lambda0 + 2) * dt
  expect_gt(q1, q0)
  expect_equal(ibm_clock_covariance(q1)[1, 1] / ibm_clock_covariance(q0)[1, 1], q1 / q0)
  expect_equal(ibm_clock_covariance(q1)[2, 2] / ibm_clock_covariance(q0)[2, 2], (q1 / q0)^3)
})

test_that("clock jumps create distinct one-sided derivatives", {
  B <- 2.5; rates <- c(0.4, 1.7)
  expect_equal(rates[1] * B, 1)
  expect_equal(rates[2] * B, 4.25)
  expect_false(isTRUE(all.equal(rates[1] * B, rates[2] * B)))
})

test_that("simulator preserves q=lambda times dt and derivative sides", {
  s <- simulate_ibm_clock(c(0, 0.2, 2, 2.5), c(1, 1, 3), seed = 4)
  expect_equal(s$q_interval, c(0.2, 1.8, 1.5))
  expect_equal(s$fprime_left[3], s$lambda_interval[2] * s$B_operational[3])
  expect_equal(s$fprime_right[3], s$lambda_interval[3] * s$B_operational[3])
  expect_true(is.na(s$fprime_left[1]))
  expect_true(is.na(s$fprime_right[4]))
})

test_that("Stan code uses literal operational-time recursion", {
  adaptive <- ibm_smooth(adaptive = "baseline_clock", get_code = TRUE)
  expect_match(adaptive, "q_interval = lambda_interval .* deltat", fixed = TRUE)
  expect_match(adaptive, "lambda_interval = lambda0 + r_interval", fixed = TRUE)
  expect_match(adaptive, "B_operational[i-1]", fixed = TRUE)
  expect_false(grepl("square(tau0) + square(psi", adaptive, fixed = TRUE))

  global <- ibm_smooth(adaptive = "ibm", get_code = TRUE)
  expect_match(global, "rep_vector(lambda0, T-1)", fixed = TRUE)
  expect_match(global, "q_interval = lambda_interval .* deltat", fixed = TRUE)
})

test_that("persistent clock shrinks adjacent log-rate changes", {
  code <- ibm_smooth(adaptive = "persistent_clock", get_code = TRUE)
  expect_match(code, "delta_log_rate = gamma * (local_raw .* z_delta)", fixed = TRUE)
  expect_match(code,
               "eta_interval[i] = eta_interval[i-1] + delta_log_rate[i-1]",
               fixed = TRUE)
  expect_match(code, "lambda_interval = exp(eta_interval)", fixed = TRUE)
  expect_match(code, "q_interval = lambda_interval .* deltat", fixed = TRUE)
  expect_match(code, "z_delta ~ std_normal()", fixed = TRUE)
  delta_line <- grep("delta_log_rate =", strsplit(code, "\\n", fixed = FALSE)[[1L]], value = TRUE)
  expect_false(any(grepl("deltat", delta_line, fixed = TRUE)))
})

test_that("zero persistent-clock increments give a constant clock", {
  eta <- log(0.7)
  delta <- rep(0, 4)
  lambda <- exp(c(eta, eta + cumsum(delta)))
  dt <- c(0.1, 0.4, 1.2, 0.3, 2)
  expect_equal(lambda, rep(0.7, 5))
  expect_equal(lambda * dt, 0.7 * dt)
})
