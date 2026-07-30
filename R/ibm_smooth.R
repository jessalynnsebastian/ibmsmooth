.ibm_stan_path <- function(adaptive,
                           parameterization = c("noncentered", "centered"),
                           adaptive_prior = c("horseshoe",
                                              "regularized_horseshoe")) {
  parameterization <- match.arg(parameterization)
  adaptive_prior <- match.arg(adaptive_prior)
  filename <- if (isTRUE(adaptive)) {
    if (parameterization == "centered" &&
        adaptive_prior == "regularized_horseshoe") {
      "ibm_adaptive_centered_regularized.stan"
    } else if (parameterization == "centered") {
      "ibm_adaptive_centered.stan"
    } else if (adaptive_prior == "regularized_horseshoe") {
      "ibm_adaptive_regularized.stan"
    } else {
      "ibm_adaptive.stan"
    }
  } else {
    if (parameterization == "centered") "ibm_centered.stan" else "ibm.stan"
  }
  installed <- system.file("stan", filename, package = "ibmsmooth")
  if (nzchar(installed)) return(installed)
  source_path <- file.path("inst", "stan", filename)
  if (file.exists(source_path)) return(source_path)
  stop("Could not locate ", filename, ".", call. = FALSE)
}

.ibm_stan_code <- function(
    adaptive, smoothing_prior,
    parameterization = c("noncentered", "centered"),
    adaptive_prior = c("horseshoe", "regularized_horseshoe")) {
  parameterization <- match.arg(parameterization)
  adaptive_prior <- match.arg(adaptive_prior)
  code <- paste(
    readLines(
      .ibm_stan_path(adaptive, parameterization, adaptive_prior),
      warn = FALSE
    ),
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

.ibm_reference_sd <- function(grid) {
  stochastic_sd <- sqrt(grid[-1L]^3 / 3)
  exp(mean(log(stochastic_sd)))
}

.ibm_global_scale <- function(reference_sd, upper, alpha,
                              prior = c("half_cauchy", "half_normal")) {
  prior <- match.arg(prior)
  quantile <- if (prior == "half_normal") {
    stats::qnorm(1 - alpha / 2)
  } else {
    tan(pi / 2 * (1 - alpha))
  }
  upper / (reference_sd * quantile)
}

.ibm_fit <- function(t, y, infer_at, adaptive, smoothing_prior, global_prior,
                     global_upper, global_alpha,
                     adaptive_prior, slab_scale, slab_df,
                     regularized_retry = "never",
                     log_sigma, initial_sd, iter, chains, cores, init = "random",
                     max_treedepth, adapt_delta, get_code, prior_only = FALSE,
                     parameterization = c("noncentered", "centered"), ...) {
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("adaptive must be TRUE or FALSE.", call. = FALSE)
  }
  parameterization <- match.arg(parameterization)
  global_prior <- match.arg(global_prior, c("half_cauchy", "half_normal"))
  adaptive_prior <- match.arg(
    adaptive_prior, c("horseshoe", "regularized_horseshoe")
  )
  regularized_retry <- match.arg(regularized_retry, c("ask", "never"))
  code <- .ibm_stan_code(
    adaptive, smoothing_prior, parameterization, adaptive_prior
  )
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
  for (value in c(global_upper, initial_sd, slab_scale, slab_df)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(paste(
        "global_upper, initial_sd, slab_scale, and slab_df must be",
        "positive finite scalars."
      ),
           call. = FALSE)
    }
  }
  if (!is.numeric(global_alpha) || length(global_alpha) != 1L ||
      !is.finite(global_alpha) || global_alpha <= 0 || global_alpha >= 1) {
    stop("global_alpha must be a finite scalar strictly between 0 and 1.",
         call. = FALSE)
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
  t_min <- min(grid_raw)
  time_scale <- diff(range(grid_raw))
  grid <- (grid_raw - t_min) / time_scale
  obs_idx <- match((t_raw - t_min) / time_scale, grid)
  reference_sd <- .ibm_reference_sd(grid)
  calibrated_global_scale <- .ibm_global_scale(
    reference_sd, global_upper, global_alpha, global_prior
  )

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
  if (isTRUE(adaptive)) {
    stan_data$global_scale <- calibrated_global_scale
    stan_data$global_prior <- match(global_prior,
                                    c("half_cauchy", "half_normal"))
  }
  if (isTRUE(adaptive) && adaptive_prior == "regularized_horseshoe") {
    stan_data$slab_scale <- slab_scale
    stan_data$slab_df <- slab_df
  }

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
      time_scale = time_scale, dt_mean = time_scale,
      y_mean = y_mean, y_sd = y_sd,
      reference_sd = reference_sd,
      infer_at = sort(unique(as.numeric(infer_at))),
      prediction_grid_raw = sort(unique(c(grid_raw, as.numeric(infer_at))))
    ),
    adaptive = adaptive,
    adaptive_prior = if (adaptive) adaptive_prior else NULL,
    parameterization = parameterization,
    smoothing_prior = if (adaptive) NULL else smoothing_prior,
    prior_only = isTRUE(prior_only),
    fit_spec = list(
      t = t_raw, y = y_raw,
      infer_at = sort(unique(as.numeric(infer_at))),
      adaptive = adaptive, parameterization = parameterization,
      smoothing_prior = smoothing_prior,
      global_prior = global_prior, global_upper = global_upper,
      global_alpha = global_alpha,
      global_scale = calibrated_global_scale,
      reference_sd = reference_sd, adaptive_prior = adaptive_prior,
      slab_scale = slab_scale, slab_df = slab_df,
      regularized_retry = regularized_retry, log_sigma = log_sigma,
      initial_sd = initial_sd, max_treedepth = max_treedepth,
      adapt_delta = adapt_delta
    )
  )
  class(fit) <- c("ibmfit", "list")

  if (isTRUE(adaptive) && adaptive_prior == "horseshoe" &&
      !isTRUE(prior_only) && regularized_retry != "never") {
    recommendation <- check_regularized_horseshoe(fit)
    if (isTRUE(recommendation$recommended)) {
      reason <- paste(recommendation$reasons, collapse = "; ")
      message(
        "The ordinary horseshoe fit shows a failure mode that a finite slab ",
        "may help: ", reason, "."
      )
      retry <- FALSE
      if (regularized_retry == "ask" && interactive()) {
        retry <- isTRUE(utils::askYesNo(paste(
          "Refit now with adaptive_prior = \"regularized_horseshoe\"",
          "and return that fit instead?"
        )))
      } else if (regularized_retry == "ask") {
        message(
          "This is a non-interactive session, so no refit was started. ",
          "Set adaptive_prior = \"regularized_horseshoe\" to retry."
        )
      }
      if (retry) {
        regularized_fit <- .ibm_fit(
          t = t, y = y, infer_at = infer_at, adaptive = adaptive,
          smoothing_prior = smoothing_prior, global_prior = global_prior,
          global_upper = global_upper, global_alpha = global_alpha,
          adaptive_prior = "regularized_horseshoe",
          slab_scale = slab_scale, slab_df = slab_df,
          regularized_retry = "never",
          log_sigma = log_sigma, initial_sd = initial_sd,
          iter = iter, chains = chains, cores = cores, init = init,
          max_treedepth = max_treedepth, adapt_delta = adapt_delta,
          get_code = FALSE, prior_only = FALSE,
          parameterization = parameterization, ...
        )
        regularized_fit$retry_info <- list(
          triggered = TRUE,
          original_adaptive_prior = "horseshoe",
          recommendation = recommendation
        )
        return(regularized_fit)
      }
      fit$regularized_horseshoe_recommendation <- recommendation
    }
  }
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
                       parameterization = c("noncentered", "centered"),
                       smoothing_prior = "tau ~ lognormal(-2, 0.5);",
                       global_prior = c("half_cauchy", "half_normal"),
                       global_upper = 1, global_alpha = 0.05,
                       adaptive_prior = c("horseshoe",
                                          "regularized_horseshoe"),
                       slab_scale = 1, slab_df = 4,
                       regularized_retry = c("ask", "never"),
                       log_sigma = list(mu = -1, sd = 1),
                       initial_sd = 5,
                       iter = 2000, chains = 4,
                       cores = getOption("mc.cores", chains),
                       init = "random",
                       max_treedepth = 12, adapt_delta = 0.9,
                       get_code = FALSE, ...) {
  ibm(
    t = t, y = y, infer_at = infer_at, adaptive = adaptive,
    parameterization = parameterization,
    smoothing_prior = smoothing_prior, global_prior = global_prior,
    global_upper = global_upper, global_alpha = global_alpha,
    adaptive_prior = adaptive_prior, slab_scale = slab_scale,
    slab_df = slab_df, regularized_retry = regularized_retry,
    log_sigma = log_sigma, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, ...
  )
}
