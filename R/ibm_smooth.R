.ibm_stan_path <- function(adaptive, fast = FALSE) {
  filename <- if (isTRUE(adaptive)) {
    if (isTRUE(fast)) "ibm_adaptive_fast.stan" else "ibm_adaptive.stan"
  } else {
    if (isTRUE(fast)) "ibm_fast.stan" else "ibm.stan"
  }
  installed <- system.file("stan", filename, package = "ibmsmooth")
  if (nzchar(installed)) return(installed)
  source_path <- file.path("inst", "stan", filename)
  if (file.exists(source_path)) return(source_path)
  stop("Could not locate ", filename, ".", call. = FALSE)
}

.ibm_stan_code <- function(adaptive, smoothing_prior, fast = FALSE) {
  code <- paste(
    readLines(.ibm_stan_path(adaptive, fast), warn = FALSE),
    collapse = "\n"
  )
  if (!isTRUE(adaptive)) {
    if (!is.character(smoothing_prior) || length(smoothing_prior) != 1L ||
        is.na(smoothing_prior) || !nzchar(trimws(smoothing_prior))) {
      stop("smoothing_prior must be one non-empty character string.", call. = FALSE)
    }
    code <- sub("// SMOOTHING_PRIOR", smoothing_prior, code, fixed = TRUE)
  }
  code
}

.ibm_model_cache <- new.env(parent = emptyenv())

.ibm_stan_model <- function(code) {
  if (!exists(code, envir = .ibm_model_cache, inherits = FALSE)) {
    model <- rstan::stan_model(model_code = code)
    assign(code, model, envir = .ibm_model_cache)
  }
  get(code, envir = .ibm_model_cache, inherits = FALSE)
}

.ibm_fit <- function(t, y, infer_at, adaptive, smoothing_prior, global_scale,
                     log_sigma, initial_sd, iter, chains, cores, init = "random",
                     max_treedepth, adapt_delta, get_code, prior_only = FALSE,
                     fast = FALSE, ...) {
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("adaptive must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("fast must be TRUE or FALSE.", call. = FALSE)
  }
  code <- .ibm_stan_code(adaptive, smoothing_prior, fast)
  if (isTRUE(get_code) || (is.null(t) && is.null(y))) return(code)

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required to fit ibmsmooth models.", call. = FALSE)
  }
  if (!is.numeric(t) || !is.numeric(y) || length(t) != length(y) ||
      length(t) < 2L || any(!is.finite(t)) || any(!is.finite(y))) {
    stop("t and y must be finite numeric vectors of equal length.", call. = FALSE)
  }
  if (!is.null(infer_at) &&
      (!is.numeric(infer_at) || any(!is.finite(infer_at)))) {
    stop("infer_at must be NULL or a finite numeric vector.", call. = FALSE)
  }
  if (!is.list(log_sigma) || !is.numeric(log_sigma$mu) ||
      length(log_sigma$mu) != 1L || !is.finite(log_sigma$mu) ||
      !is.numeric(log_sigma$sd) || length(log_sigma$sd) != 1L ||
      !is.finite(log_sigma$sd) || log_sigma$sd <= 0) {
    stop("log_sigma must contain a finite mu and a positive finite sd.", call. = FALSE)
  }
  for (value in c(global_scale, initial_sd)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("global_scale and initial_sd must be positive finite scalars.",
           call. = FALSE)
    }
  }

  t_raw <- as.numeric(t)
  y_raw <- as.numeric(y)
  # Only observed locations belong in the Stan state vector. States requested
  # through infer_at are recovered afterward from exact IBM bridges.
  grid_raw <- sort(unique(t_raw))
  if (length(grid_raw) < 2L) {
    stop("At least two unique locations are required.", call. = FALSE)
  }
  y_mean <- mean(y_raw)
  y_sd <- stats::sd(y_raw)
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("y must have positive finite standard deviation.", call. = FALSE)
  }
  dt_mean <- mean(diff(grid_raw))
  t_min <- min(grid_raw)
  grid <- (grid_raw - t_min) / dt_mean
  obs_idx <- match((t_raw - t_min) / dt_mean, grid)

  stan_data <- list(
    N_obs = length(y_raw),
    T = length(grid),
    obs_time_idx = obs_idx,
    deltat = diff(grid),
    y_obs = (y_raw - y_mean) / y_sd,
    log_sigma_mu = log_sigma$mu,
    log_sigma_sd = log_sigma$sd,
    initial_sd = initial_sd,
    prior_only = as.integer(isTRUE(prior_only))
  )
  if (isTRUE(adaptive)) stan_data$global_scale <- global_scale

  stan_model <- .ibm_stan_model(code)
  stanfit <- rstan::sampling(
    object = stan_model, data = stan_data, iter = iter, chains = chains,
    cores = cores, init = init,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta),
    ...
  )

  fit <- list(
    stanfit = stanfit,
    data = list(
      t_raw = t_raw, y_raw = y_raw, time_grid_raw = grid_raw,
      time_grid = grid, obs_time_idx = obs_idx, t_min = t_min,
      dt_mean = dt_mean, y_mean = y_mean, y_sd = y_sd,
      infer_at = sort(unique(as.numeric(infer_at))),
      prediction_grid_raw = sort(unique(c(grid_raw, as.numeric(infer_at))))
    ),
    adaptive = adaptive,
    fast = isTRUE(fast),
    smoothing_prior = if (adaptive) NULL else smoothing_prior,
    prior_only = isTRUE(prior_only),
    fit_spec = list(
      t = t_raw, y = y_raw,
      infer_at = sort(unique(as.numeric(infer_at))),
      adaptive = adaptive, fast = isTRUE(fast),
      smoothing_prior = smoothing_prior,
      global_scale = global_scale, log_sigma = log_sigma,
      initial_sd = initial_sd, max_treedepth = max_treedepth,
      adapt_delta = adapt_delta
    )
  )
  class(fit) <- c("ibmfit", "list")
  fit
}

#' Fit an IBM smoother using Stan
#'
#' `ibm_smooth()` is an alias for [ibm()] retained as the package's
#' lower-level fitting name.
#'
#' @inheritParams ibm
#' @return An object of class `ibmfit`, or Stan code when `get_code = TRUE`.
#' @export
ibm_smooth <- function(t = NULL, y = NULL, infer_at = NULL, adaptive = FALSE,
                       fast = FALSE,
                       smoothing_prior = "tau ~ lognormal(-2, 0.5);",
                       global_scale = 0.1,
                       log_sigma = list(mu = -1, sd = 1),
                       initial_sd = 5,
                       iter = 2000, chains = 4,
                       cores = getOption("mc.cores", chains),
                       init = "random",
                       max_treedepth = 12, adapt_delta = 0.9,
                       get_code = FALSE, ...) {
  ibm(
    t = t, y = y, infer_at = infer_at, adaptive = adaptive, fast = fast,
    smoothing_prior = smoothing_prior, global_scale = global_scale,
    log_sigma = log_sigma, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, ...
  )
}
