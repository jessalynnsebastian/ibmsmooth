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

.ibm_families <- c(
  "gaussian", "student_t", "bernoulli", "binomial", "poisson",
  "negative_binomial", "lognormal", "gamma", "beta", "exponential"
)

.ibm_likelihood_code <- function(family) {
  family <- match.arg(family, .ibm_families)
  scale_data <- paste(
    "vector[N_obs] y_obs;", "real log_sigma_mu;",
    "real<lower=0> log_sigma_sd;", sep = "\n  "
  )
  scale_parts <- list(
    data = scale_data,
    parameters = "real log_sigma_raw;",
    transformed = paste0(
      "real<lower=0> sigma = exp(log_sigma_mu + ",
      "log_sigma_sd * log_sigma_raw);"
    ),
    priors = "log_sigma_raw ~ std_normal();"
  )
  phi_parts <- list(
    data = paste(
      "real log_phi_mu;", "real<lower=0> log_phi_sd;", sep = "\n  "
    ),
    parameters = "real log_phi_raw;",
    transformed = paste0(
      "real<lower=0> phi = exp(log_phi_mu + log_phi_sd * log_phi_raw);"
    ),
    priors = "log_phi_raw ~ std_normal();"
  )
  empty <- list(data = "", parameters = "", transformed = "", priors = "")
  parts <- switch(
    family,
    gaussian = scale_parts,
    student_t = {
      scale_parts$data <- paste(scale_parts$data,
                                "real<lower=0> student_df;", sep = "\n  ")
      scale_parts
    },
    lognormal = scale_parts,
    gamma = phi_parts,
    beta = phi_parts,
    negative_binomial = phi_parts,
    empty
  )
  extra_data <- switch(
    family,
    bernoulli = "array[N_obs] int<lower=0, upper=1> y_int;",
    binomial = paste(
      "array[N_obs] int<lower=0> y_int;",
      "array[N_obs] int<lower=0> trials;", sep = "\n  "
    ),
    poisson = paste(
      "array[N_obs] int<lower=0> y_int;",
      "vector<lower=0>[N_obs] exposure;", sep = "\n  "
    ),
    negative_binomial = paste(
      "array[N_obs] int<lower=0> y_int;",
      "vector<lower=0>[N_obs] exposure;", sep = "\n  "
    ),
    exponential = "vector<lower=0>[N_obs] y_obs;",
    gamma = "vector<lower=0>[N_obs] y_obs;",
    beta = "vector<lower=0, upper=1>[N_obs] y_obs;",
    ""
  )
  if (nzchar(extra_data)) {
    parts$data <- paste(c(parts$data, extra_data)[nzchar(c(parts$data,
                                                           extra_data))],
                        collapse = "\n  ")
  }
  likelihood <- switch(
    family,
    gaussian = "y_obs ~ normal(f[obs_time_idx], sigma);",
    student_t = "y_obs ~ student_t(student_df, f[obs_time_idx], sigma);",
    bernoulli = "y_int ~ bernoulli_logit(f[obs_time_idx]);",
    binomial = "y_int ~ binomial_logit(trials, f[obs_time_idx]);",
    poisson = paste0(
      "y_int ~ poisson_log(f[obs_time_idx] + log(exposure));"
    ),
    negative_binomial = paste0(
      "y_int ~ neg_binomial_2_log(f[obs_time_idx] + log(exposure), phi);"
    ),
    lognormal = "y_obs ~ lognormal(f[obs_time_idx], sigma);",
    gamma = paste0(
      "for (n in 1:N_obs) y_obs[n] ~ gamma(phi, ",
      "phi / exp(f[obs_time_idx[n]]));"
    ),
    beta = paste0(
      "for (n in 1:N_obs) { real mu = inv_logit(f[obs_time_idx[n]]); ",
      "y_obs[n] ~ beta(mu * phi, (1 - mu) * phi); }"
    ),
    exponential = "y_obs ~ exponential(exp(-f[obs_time_idx]));"
  )
  c(parts, list(likelihood = likelihood))
}

.ibm_stan_code <- function(
    adaptive, smoothing_prior,
    parameterization = c("noncentered", "centered"),
    adaptive_prior = c("horseshoe", "regularized_horseshoe"),
    family = .ibm_families) {
  parameterization <- match.arg(parameterization)
  adaptive_prior <- match.arg(adaptive_prior)
  family <- match.arg(family, .ibm_families)
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
  observation <- .ibm_likelihood_code(family)
  replacements <- c(
    OBSERVATION_DATA = "data", OBSERVATION_PARAMETERS = "parameters",
    OBSERVATION_TRANSFORMED_PARAMETERS = "transformed",
    OBSERVATION_PRIORS = "priors", OBSERVATION_LIKELIHOOD = "likelihood"
  )
  for (marker in names(replacements)) {
    code <- sub(paste0("// ", marker),
                observation[[replacements[[marker]]]], code, fixed = TRUE)
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

.ibm_reference_sd <- function(grid = NULL) {
  # On the standardized domain, sd{g(t)} = t^(3/2) / sqrt(3), so its
  # grid-independent maximum over [0, 1] is attained at t = 1.
  1 / sqrt(3)
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
                     family, trials, exposure, student_df, log_sigma, log_phi,
                     initial_sd, iter, chains, cores, init = "random",
                     max_treedepth, adapt_delta, get_code, prior_only = FALSE,
                     parameterization = c("noncentered", "centered"),
                     print_code = FALSE, stan_file = NULL, ...) {
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) {
    stop("adaptive must be TRUE or FALSE.", call. = FALSE)
  }
  parameterization <- match.arg(parameterization)
  family <- match.arg(family, .ibm_families)
  global_prior <- match.arg(global_prior, c("half_cauchy", "half_normal"))
  adaptive_prior <- match.arg(
    adaptive_prior, c("horseshoe", "regularized_horseshoe")
  )
  regularized_retry <- match.arg(regularized_retry, c("ask", "never"))
  code <- .ibm_stan_code(
    adaptive, smoothing_prior, parameterization, adaptive_prior, family
  )
  if (!is.logical(print_code) || length(print_code) != 1L ||
      is.na(print_code)) {
    stop("print_code must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(stan_file)) {
    if (!is.character(stan_file) || length(stan_file) != 1L ||
        is.na(stan_file) || !nzchar(stan_file)) {
      stop("stan_file must be NULL or one non-empty file path.", call. = FALSE)
    }
    writeLines(code, stan_file, useBytes = TRUE)
  }
  if (isTRUE(print_code)) cat(code, "\n")
  if (isTRUE(get_code)) return(code)
  if (is.null(t) && is.null(y)) {
    if (isTRUE(print_code) || !is.null(stan_file)) return(invisible(code))
    return(code)
  }

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
  if (family %in% c("gaussian", "student_t", "lognormal") &&
      (!is.list(log_sigma) || !is.numeric(log_sigma$mu) ||
      length(log_sigma$mu) != 1L || !is.finite(log_sigma$mu) ||
      !is.numeric(log_sigma$sd) || length(log_sigma$sd) != 1L ||
      !is.finite(log_sigma$sd) || log_sigma$sd <= 0)) {
    stop("log_sigma must contain a finite mu and a positive finite sd.", call. = FALSE)
  }
  if (family %in% c("negative_binomial", "gamma", "beta") &&
      (!is.list(log_phi) || !is.numeric(log_phi$mu) ||
      length(log_phi$mu) != 1L || !is.finite(log_phi$mu) ||
      !is.numeric(log_phi$sd) || length(log_phi$sd) != 1L ||
      !is.finite(log_phi$sd) || log_phi$sd <= 0)) {
    stop("log_phi must contain a finite mu and a positive finite sd.",
         call. = FALSE)
  }
  if (family == "student_t" &&
      (!is.numeric(student_df) || length(student_df) != 1L ||
       !is.finite(student_df) || student_df <= 0)) {
    stop("student_df must be a positive finite scalar.", call. = FALSE)
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
  standardize_response <- family %in% c("gaussian", "student_t")
  if (standardize_response && (!is.finite(y_sd) || y_sd <= 0)) {
    stop("y must have positive finite standard deviation.", call. = FALSE)
  }
  if (!standardize_response) {
    y_mean <- 0
    y_sd <- 1
  }
  integer_family <- family %in% c(
    "bernoulli", "binomial", "poisson", "negative_binomial"
  )
  if (integer_family && any(y_raw < 0 | y_raw != floor(y_raw))) {
    stop("y must contain nonnegative integers for this family.",
         call. = FALSE)
  }
  if (family == "bernoulli" && any(!y_raw %in% 0:1)) {
    stop("Bernoulli observations must be zero or one.", call. = FALSE)
  }
  if (family == "binomial") {
    if (is.null(trials) || !is.numeric(trials) ||
        length(trials) != length(y_raw) || any(!is.finite(trials)) ||
        any(trials < 0 | trials != floor(trials)) || any(y_raw > trials)) {
      stop(paste(
        "trials must contain one nonnegative integer per observation and",
        "must be at least as large as y."
      ), call. = FALSE)
    }
  }
  if (family %in% c("poisson", "negative_binomial")) {
    if (!is.numeric(exposure) || !length(exposure) ||
        !(length(exposure) %in% c(1L, length(y_raw))) ||
        any(!is.finite(exposure)) || any(exposure <= 0)) {
      stop("exposure must be a positive scalar or one value per observation.",
           call. = FALSE)
    }
    exposure <- rep(as.numeric(exposure), length.out = length(y_raw))
  }
  if (family %in% c("lognormal", "gamma", "exponential") &&
      any(y_raw <= 0)) {
    stop("y must be strictly positive for this family.", call. = FALSE)
  }
  if (family == "beta" && any(y_raw <= 0 | y_raw >= 1)) {
    stop("Beta observations must lie strictly between zero and one.",
         call. = FALSE)
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
    initial_sd = initial_sd,
    prior_only = as.integer(isTRUE(prior_only))
  )
  if (family %in% c("gaussian", "student_t")) {
    stan_data$y_obs <- (y_raw - y_mean) / y_sd
  } else if (family %in% c("lognormal", "gamma", "beta", "exponential")) {
    stan_data$y_obs <- y_raw
  } else {
    stan_data$y_int <- as.integer(y_raw)
  }
  if (family %in% c("gaussian", "student_t", "lognormal")) {
    stan_data$log_sigma_mu <- log_sigma$mu
    stan_data$log_sigma_sd <- log_sigma$sd
  }
  if (family %in% c("negative_binomial", "gamma", "beta")) {
    stan_data$log_phi_mu <- log_phi$mu
    stan_data$log_phi_sd <- log_phi$sd
  }
  if (family == "student_t") stan_data$student_df <- student_df
  if (family == "binomial") stan_data$trials <- as.integer(trials)
  if (family %in% c("poisson", "negative_binomial")) {
    stan_data$exposure <- exposure
  }
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
    family = family,
    adaptive_prior = if (adaptive) adaptive_prior else NULL,
    parameterization = parameterization,
    smoothing_prior = if (adaptive) NULL else smoothing_prior,
    prior_only = isTRUE(prior_only),
    fit_spec = list(
      t = t_raw, y = y_raw,
      infer_at = sort(unique(as.numeric(infer_at))),
      adaptive = adaptive, family = family, trials = trials,
      exposure = exposure, student_df = student_df,
      parameterization = parameterization,
      smoothing_prior = smoothing_prior,
      global_prior = global_prior, global_upper = global_upper,
      global_alpha = global_alpha,
      global_scale = calibrated_global_scale,
      reference_sd = reference_sd, adaptive_prior = adaptive_prior,
      slab_scale = slab_scale, slab_df = slab_df,
      regularized_retry = regularized_retry, log_sigma = log_sigma,
      log_phi = log_phi,
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
          family = family, trials = trials, exposure = exposure,
          student_df = student_df,
          smoothing_prior = smoothing_prior, global_prior = global_prior,
          global_upper = global_upper, global_alpha = global_alpha,
          adaptive_prior = "regularized_horseshoe",
          slab_scale = slab_scale, slab_df = slab_df,
          regularized_retry = "never",
          log_sigma = log_sigma, log_phi = log_phi,
          initial_sd = initial_sd,
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
                       family = c("gaussian", "student_t", "bernoulli",
                                  "binomial", "poisson", "negative_binomial",
                                  "lognormal", "gamma", "beta",
                                  "exponential"),
                       trials = NULL, exposure = 1, student_df = 4,
                       parameterization = c("noncentered", "centered"),
                       smoothing_prior = "tau ~ lognormal(-2, 0.5);",
                       global_prior = c("half_cauchy", "half_normal"),
                       global_upper = 1, global_alpha = 0.05,
                       adaptive_prior = c("horseshoe",
                                          "regularized_horseshoe"),
                       slab_scale = 1, slab_df = 4,
                       regularized_retry = c("ask", "never"),
                       log_sigma = list(mu = -1, sd = 1),
                       log_phi = list(mu = 0, sd = 1),
                       initial_sd = 5,
                       iter = 2000, chains = 4,
                       cores = getOption("mc.cores", chains),
                       init = "random",
                       max_treedepth = 12, adapt_delta = 0.9,
                       get_code = FALSE, print_code = FALSE,
                       stan_file = NULL, ...) {
  ibm(
    t = t, y = y, infer_at = infer_at, adaptive = adaptive,
    family = family, trials = trials, exposure = exposure,
    student_df = student_df,
    parameterization = parameterization,
    smoothing_prior = smoothing_prior, global_prior = global_prior,
    global_upper = global_upper, global_alpha = global_alpha,
    adaptive_prior = adaptive_prior, slab_scale = slab_scale,
    slab_df = slab_df, regularized_retry = regularized_retry,
    log_sigma = log_sigma, log_phi = log_phi, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, print_code = print_code, stan_file = stan_file, ...
  )
}
