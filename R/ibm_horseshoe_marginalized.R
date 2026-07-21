#' Fit a marginalized Gaussian horseshoe IBM
#'
#' Fits the same Gaussian-data model as the current noncentered horseshoe Stan
#' model, but integrates the two-dimensional IBM state out with a Kalman filter
#' during HMC. Exact conditional state trajectories are reconstructed afterward
#' with a simulation smoother.
#'
#' The model is also available through [ibm()] by selecting the marginalized
#' horseshoe engine.
#'
#' @param t Numeric vector of observation times.
#' @param y Numeric vector of observations.
#' @param infer_at Optional additional times at which state draws are required.
#' @param log_sigma Prior specification used by the existing Stan backend.
#' @param zeta Global horseshoe scale.
#' @param initial_sd Scale of the correlated unit-time IBM initial-state prior.
#' @param n_state_draws Number of conditional state trajectories to reconstruct.
#'   Set to zero to retain only hyperparameter draws.
#' @param seed Optional random seed used by Stan and the simulation smoother.
#' @param iter,chains,cores Stan sampling controls.
#' @param max_treedepth,adapt_delta Stan HMC controls.
#' @param ... Additional arguments passed to [rstan::stan()].
#'
#' @return An object of class `ibmfit` containing the marginalized Stan fit and
#'   cached posterior state draws.
#' @export
ibm_horseshoe_marginalized <- function(
    t, y, infer_at = NULL,
    log_sigma = list(mu = -1, sd = 1), zeta = 0.1, initial_sd = 5,
    n_state_draws = 1000L, seed = NULL,
    iter = 2000, chains = 4,
    cores = getOption("mc.cores", chains),
    max_treedepth = 12, adapt_delta = 0.9, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required for the marginalized horseshoe model.")
  }
  n_state_draws <- as.integer(n_state_draws)
  if (length(n_state_draws) != 1L || is.na(n_state_draws) || n_state_draws < 0L)
    stop("n_state_draws must be a non-negative integer")
  prep <- .ibmsmooth_marginalized_data(t, y, infer_at, log_sigma, zeta, initial_sd)
  stan_file <- system.file(
    "stan", "ibm_adaptive_horseshoe_gaussian_marginalized.stan",
    package = "ibmsmooth"
  )
  if (!nzchar(stan_file)) {
    candidate <- file.path("inst", "stan", "ibm_adaptive_horseshoe_gaussian_marginalized.stan")
    if (file.exists(candidate)) stan_file <- candidate
  }
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop("Could not find the marginalized horseshoe Stan model.")
  }

  stan_args <- list(
    file = stan_file, data = prep$stan_data, iter = iter, chains = chains,
    cores = cores,
    control = list(max_treedepth = max_treedepth, adapt_delta = adapt_delta)
  )
  if (!is.null(seed)) stan_args$seed <- seed
  stanfit <- do.call(rstan::stan, c(stan_args, list(...)))

  fit <- list(
    stanfit = stanfit,
    data = prep$data,
    marginalized = TRUE,
    stan_model = "horseshoe_gaussian_marginalized",
    state_samples_scaled = NULL
  )
  class(fit) <- c("ibmfit", class(fit))

  if (n_state_draws > 0L) {
    hp <- rstan::extract(stanfit, pars = c("sigma", "tau"), permuted = TRUE)
    n_available <- length(hp$sigma)
    keep <- if (n_state_draws >= n_available) seq_len(n_available) else
      unique(round(seq(1, n_available, length.out = n_state_draws)))
    if (!is.null(seed)) set.seed(seed + 1L)
    draws <- lapply(keep, function(k) {
      .ibmsmooth_ffbs(
        prep$stan_data$y_obs, prep$stan_data$obs_time_idx,
        prep$stan_data$deltat, hp$sigma[k], hp$tau[k, ], initial_sd
      )
    })
    fit$state_samples_scaled <- list(
      f = do.call(rbind, lapply(draws, function(x) x[2L, ])),
      fprime = do.call(rbind, lapply(draws, function(x) x[1L, ])),
      hyperparameter_draw = keep
    )
  }
  fit
}

.ibmsmooth_marginalized_data <- function(t, y, infer_at, log_sigma, zeta,
                                         initial_sd = 5) {
  if (!is.numeric(t) || !is.numeric(y) || length(t) != length(y))
    stop("t and y must be numeric vectors of the same length")
  if (length(t) < 2L || any(!is.finite(t)) || any(!is.finite(y)))
    stop("t and y must contain at least two finite observations")
  if (is.null(log_sigma$mu) || is.null(log_sigma$sd) || log_sigma$sd <= 0)
    stop("log_sigma must contain mu and a positive sd")
  if (!is.numeric(zeta) || length(zeta) != 1L || !is.finite(zeta) || zeta <= 0)
    stop("zeta must be a positive finite scalar")
  if (!is.numeric(initial_sd) || length(initial_sd) != 1L ||
      !is.finite(initial_sd) || initial_sd <= 0)
    stop("initial_sd must be a positive finite scalar")

  time_grid_raw <- sort(unique(c(as.numeric(t), as.numeric(infer_at))))
  if (length(time_grid_raw) < 2L) stop("At least two unique grid times are required")
  y_mean <- mean(y)
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd <= 0) stop("y must have positive finite standard deviation")
  dt_mean <- mean(diff(time_grid_raw))
  t_min <- min(time_grid_raw)
  time_grid <- (time_grid_raw - t_min) / dt_mean
  obs_time_idx <- match((t - t_min) / dt_mean, time_grid)
  ord <- order(obs_time_idx)
  dt <- diff(time_grid)

  list(
    stan_data = list(
      N_obs = length(y), T = length(time_grid),
      obs_time_idx = obs_time_idx[ord], deltat = dt,
      y_obs = ((y - y_mean) / y_sd)[ord],
      log_sigma_mu = log_sigma$mu, log_sigma_sd = log_sigma$sd,
      zeta = zeta, initial_sd = initial_sd
    ),
    data = list(
      y_raw = as.numeric(y), y_mean = y_mean, y_sd = y_sd,
      t_raw = as.numeric(t), t_min = t_min, dt_mean = dt_mean,
      time_grid_raw = time_grid_raw, time_grid = time_grid,
      obs_time_idx = obs_time_idx,
      initial_sd = initial_sd,
      regular = all(abs(dt - dt[1L]) < sqrt(.Machine$double.eps))
    )
  )
}

.ibmsmooth_psd_draw <- function(mean, cov) {
  cov <- (cov + t(cov)) / 2
  ev <- eigen(cov, symmetric = TRUE)
  values <- pmax(ev$values, 0)
  as.numeric(mean + ev$vectors %*% (sqrt(values) * stats::rnorm(length(mean))))
}

.ibmsmooth_ffbs <- function(y, obs_time_idx, dt, sigma, tau, initial_sd = 5) {
  T <- length(dt) + 1L
  mf <- vector("list", T); Cf <- vector("list", T)
  ap <- vector("list", T); Rp <- vector("list", T)
  m <- c(0, 0)
  P <- initial_sd^2 * matrix(c(1, 0.5, 0.5, 1 / 3), 2L, 2L)

  for (i in seq_len(T)) {
    if (i > 1L) {
      h <- dt[i - 1L]
      F <- matrix(c(1, h, 0, 1), 2L, 2L)
      S <- matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2L, 2L)
      m <- as.numeric(F %*% m)
      P <- F %*% P %*% t(F) + tau[i - 1L]^2 * S
    }
    ap[[i]] <- m; Rp[[i]] <- (P + t(P)) / 2
    for (n in which(obs_time_idx == i)) {
      v <- y[n] - m[2L]
      q <- P[2L, 2L] + sigma^2
      K <- P[, 2L] / q
      m <- m + K * v
      P <- P - tcrossprod(K, P[2L, ])
      P <- (P + t(P)) / 2
    }
    mf[[i]] <- m; Cf[[i]] <- P
  }

  state <- matrix(NA_real_, 2L, T)
  state[, T] <- .ibmsmooth_psd_draw(mf[[T]], Cf[[T]])
  if (T > 1L) for (i in (T - 1L):1L) {
    h <- dt[i]
    F <- matrix(c(1, h, 0, 1), 2L, 2L)
    J <- Cf[[i]] %*% t(F) %*% solve(Rp[[i + 1L]])
    cm <- mf[[i]] + as.numeric(J %*% (state[, i + 1L] - ap[[i + 1L]]))
    CC <- Cf[[i]] - J %*% Rp[[i + 1L]] %*% t(J)
    state[, i] <- .ibmsmooth_psd_draw(cm, CC)
  }
  state
}

#' Validate the marginalized Kalman likelihood against a dense Gaussian result
#'
#' Intended as a small-data numerical check for the marginalized model.
#'
#' @param t,y Observation times and scaled observations. Times must be unique.
#' @param sigma Observation standard deviation.
#' @param tau Vector of transition scales of length `length(t) - 1`.
#' @param initial_sd Scale of the correlated unit-time IBM initial-state prior.
#' @return The Kalman and dense log likelihoods and their difference.
#' @export
validate_marginalized_horseshoe <- function(t, y, sigma, tau, initial_sd = 5) {
  if (length(unique(t)) != length(t)) stop("Validation helper requires unique times")
  ord <- order(t); t <- t[ord]; y <- y[ord]
  dt <- diff(t) / mean(diff(t))
  T <- length(t)
  if (length(tau) != T - 1L) stop("tau must have length length(t) - 1")

  initial_cov <- initial_sd^2 * matrix(c(1, 0.5, 0.5, 1 / 3), 2L, 2L)
  m <- c(0, 0); P <- initial_cov; kalman <- 0
  for (i in seq_len(T)) {
    if (i > 1L) {
      h <- dt[i - 1L]; F <- matrix(c(1, h, 0, 1), 2L, 2L)
      S <- matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2L, 2L)
      m <- as.numeric(F %*% m); P <- F %*% P %*% t(F) + tau[i - 1L]^2 * S
    }
    q <- P[2L, 2L] + sigma^2; v <- y[i] - m[2L]
    kalman <- kalman + stats::dnorm(v, 0, sqrt(q), log = TRUE)
    K <- P[, 2L] / q; m <- m + K * v; P <- P - tcrossprod(K, P[2L, ])
  }

  G <- matrix(0, 2L * T, 2L * T)
  G[1:2, 1:2] <- initial_cov
  for (i in 2:T) {
    h <- dt[i - 1L]; F <- matrix(c(1, h, 0, 1), 2L, 2L)
    S <- matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2L, 2L) * tau[i - 1L]^2
    prev <- (2L * i - 3L):(2L * i - 2L); now <- (2L * i - 1L):(2L * i)
    G[now, seq_len(2L * (i - 1L))] <- F %*% G[prev, seq_len(2L * (i - 1L))]
    G[seq_len(2L * (i - 1L)), now] <- t(G[now, seq_len(2L * (i - 1L))])
    G[now, now] <- F %*% G[prev, prev] %*% t(F) + S
  }
  level <- seq(2L, 2L * T, by = 2L)
  V <- G[level, level] + diag(sigma^2, T)
  L <- chol(V)
  dense <- -0.5 * (T * log(2 * pi) + 2 * sum(log(diag(L))) +
                     sum(backsolve(L, y, transpose = TRUE)^2))
  c(kalman = kalman, dense = dense, difference = kalman - dense)
}
