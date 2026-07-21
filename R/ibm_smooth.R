#' Fit an integrated Brownian motion (IBM) smoother using Stan
#'
#' Fit a Bayesian integrated Brownian motion smoothing model to observations
#' (t, y) using Stan. The function rescales time and response for more
#' efficient sampling, constructs the data list for Stan, and runs
#' rstan::stan() on a centered or non-centered Stan implementation.
#'
#' This function supports several adaptive smoothing options:
#' - "nonadaptive" (default): single global smoothing parameter (tau).
#' - "rw": random-walk prior on log(tau) across segments.
#' - "horseshoe": horseshoe prior on local IBM process scales.
#' - "baseline_horseshoe": horseshoe excess roughness above a positive IBM baseline.
#' - "persistent_clock": horseshoe shrinkage on adjacent log clock-rate changes.
#' - "rhs": regularized horseshoe prior on local IBM process scales.
#' - "bridge": bridge prior on whitened IBM transition innovations (requires centered parameterization).
#'
#' @param t Numeric vector of time points.
#' @param y Numeric vector of observations at times \code{t}. Must be same
#'   length as \code{t}.
#' @param infer_at Optional numeric vector of time points to do function inference at
#'  (in addition to just the observed time points). If not \code{NULL} these will be
#'  included in the Stan data and the fitted function will be returned at these points as well.
#' @param get_code Logical; if \code{TRUE} return the Stan model code as a
#'   single string (after potentially writing to \code{write_to_file}). If
#'   \code{get_code = TRUE} and both \code{t} and \code{y} are \code{NULL} the
#'   function will return the Stan code without attempting to fit a model.
#' @param write_to_file Optional character path. If provided and
#'   \code{get_code = TRUE} the Stan code is written to this file.
#' @param adaptive Character; specifies if and how to use adaptive smoothing.
#'   Options are \code{"nonadaptive"} (default), \code{"rw"}, \code{"horseshoe"},
#'   \code{"baseline_clock"}, \code{"persistent_clock"},
#'   \code{"curvature_horseshoe"}, \code{"rw"}, \code{"rhs"}, or
#'   \code{"bridge"}. Older model names remain available as aliases.
#' @param log_sigma List with elements \code{mu} and \code{sd} giving the prior
#'   mean and standard deviation for \code{log(sigma)}. Defaults to
#'   \code{list(mu = -1, sd = 1)}.
#' @param log_tau List with elements \code{mu} and \code{sd} giving the prior
#'   mean and standard deviation for \code{log(tau)}. Used for
#'   \code{adaptive = "nonadaptive"}, \code{"rw"},
#'   \code{"baseline_clock"}, and \code{"persistent_clock"}. For clock
#'   models this controls the prior for the baseline or initial clock rate.
#'   Defaults to
#'   \code{list(mu = -2, sd = 0.5)}.
#' @param zeta Numeric; global shrinkage scale used by horseshoe,
#'   baseline-clock, persistent-clock, rhs and bridge priors. For the
#'   persistent clock it controls the global scale of log-rate changes.
#'   Default \code{0.1}.
#' @param slab_scale Numeric; slab scale for the regularized horseshoe (rhs).
#'   Default \code{2}.
#' @param alpha Numeric in (0,2]; exponent for bridge prior (only used when
#'   \code{adaptive = "bridge"}). Default \code{0.5}.
#' @param initial_sd Prior scale for the correlated initial state. The Stan
#'   backend now uses the same unit-time IBM initial covariance as the INLA
#'   backend: for \code{c(f_1, fprime_1)}, the covariance is
#'   \code{initial_sd^2 * matrix(c(1/3, 1/2, 1/2, 1), 2, 2)} on the
#'   scaled response/time scale.
#' @param centered Logical; if \code{TRUE} use the centered Stan
#'   parameterization. If \code{NULL} a sensible default is chosen:
#'   centered is \code{TRUE} when \code{adaptive == "nonadaptive"} (except
#'   that \code{"bridge"} forces centered = TRUE). You can alternatively pass
#'   \code{noncentered} (logical) to explicitly request the non-centered form.
#' @param noncentered Deprecated alias for specifying non-centered parameterization.
#' @param iter Integer; total number of iterations for each chain (including
#'   warmup). Default \code{2000}.
#' @param chains Integer; number of MCMC chains. Default \code{4}.
#' @param cores Integer; number of parallel processes to use for Stan chains. Defaults to \code{getOption("mc.cores", chains)}.
#' @param max_treedepth Integer; maximum tree depth for Stan's NUTS sampler.
#'   Default \code{12}.
#' @param adapt_delta Numeric in (0,1); target average acceptance probability
#'   for Stan's NUTS. Default \code{0.85}.
#' @param ... Additional arguments passed on to \code{rstan::stan()}.
#'
#' @details
#' The function rescales \code{y} to zero mean and unit variance and rescales
#' time so that the mean time step is near 1. A flag \code{regular} is
#' computed internally to indicate whether time points are equally spaced. The
#' transition covariance for interval \code{i} is \code{tau_i^2 S_i}, where
#' \code{S_i} is the IBM covariance determined by \code{deltat[i]}. INLA uses
#' the equivalent precision parameterization \code{lambda_i = 1 / tau_i^2}.
#'
#' If \code{get_code = TRUE} the Stan code used for the chosen combination of
#' \code{adaptive} and \code{centered} is returned as a single string; if
#' \code{write_to_file} is provided the code will also be written to that path.
#'
#' @return An object of class \code{ibmfit} (a list) with components:
#'   \item{stanfit}{The \code{stanfit} object returned by \code{rstan::stan()}.}
#'   \item{data}{A list containing the original data and scaling constants
#'     used: \code{y_raw}, \code{y_mean}, \code{y_sd}, \code{t_raw},
#'     \code{t_min}, and \code{dt_mean}.}
#'
#' @seealso \code{\link[rstan]{stan}}
#'
#' @export
ibm_smooth <- function(t = NULL, y = NULL,
       infer_at = NULL, # if not null, a vector of time points to do function inference at (in addition to just the observed time points)
       # if t and y are null, we're going to be printing the code
       # or writing it to a file
       get_code = FALSE,
       write_to_file = NULL,
       # specify adaptivity
       adaptive = "nonadaptive", # constant or baseline stochastic clock
       # priors for hyperparameters
       # for nonadaptive and rw
       log_sigma = list(mu = -1, sd = 1),
       log_tau = list(mu = -2, sd = 0.5),
       # for horseshoe and rhs
       zeta = 0.1,
       slab_scale = 2,
       # for bridge
       alpha = 0.5,
       initial_sd = 5,
       centered = NULL,
       noncentered = NULL,
       iter = 2000, chains = 4,
       cores = getOption("mc.cores", chains),
       max_treedepth = 12, adapt_delta = 0.9, ...) {
  adaptive <- match.arg(adaptive, c("nonadaptive", "ibm", "baseline_clock",
                                    "adaptive_clock", "baseline_horseshoe",
                                    "persistent_clock",
                                    "curvature_horseshoe", "rw", "horseshoe",
                                    "rhs", "bridge"))
  if (adaptive == "ibm") adaptive <- "nonadaptive"
  if (adaptive %in% c("baseline_clock", "adaptive_clock")) adaptive <- "baseline_horseshoe"
  if (adaptive == "curvature_horseshoe") adaptive <- "horseshoe"
  # determine centered / noncentered choice
  if (is.null(centered) && is.null(noncentered)) {
    # if neither given, default to centered
    centered <- FALSE
  } else if (!is.null(centered) && !is.null(noncentered)) {
    warning("`noncentered` and  `centered` both provided; using `centered`.")
    centered <- as.logical(centered) # prefer centered if both given
  } else if (!is.null(centered)) {
    centered <- as.logical(centered)
  } else {
    centered <- !as.logical(noncentered)
  }
  if (!is.logical(centered) || length(centered) != 1 || is.na(centered)) {
    stop("`centered` or `noncentered` must be TRUE or FALSE")
  }
  if (adaptive == "bridge") {
    centered <- TRUE # bridge must be centered
  }
  if (adaptive %in% c("nonadaptive", "baseline_horseshoe", "persistent_clock") && centered) {
    warning("Literal stochastic-clock models use the non-centered exact transition; ignoring `centered = TRUE`.")
    centered <- FALSE
  }
  # choose appropriate stan file
  file <- "ibm"
  if (adaptive != "nonadaptive") {
    file <- paste0(file, "_adaptive_", adaptive)
  } else {
    file <- paste0(file, "_smooth")
  }
  if (centered) {
    file <- paste0(file, "_centered.stan")
  } else {
    file <- paste0(file, ".stan")
  }
  # if get_code is TRUE, return the stan code as a string
  if (is.null(t) && is.null(y)) get_code <- TRUE
  # or write to file if specified
  if (get_code) {
    stan_code <- readLines(system.file("stan", file, package = "ibmsmooth"))
    stan_code <- paste(stan_code, collapse = "\n")
    if (!is.null(write_to_file)) {
      writeLines(stan_code, con = write_to_file)
    }
    # if t and y are NULL, just return the code
    if (is.null(t) && is.null(y)) {
      return(stan_code)
    }
  }

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required for method = 'stan'. Install rstan or use method = 'inla'.")
  }

  # check inputs
  if (!is.numeric(initial_sd) || length(initial_sd) != 1 ||
      !is.finite(initial_sd) || initial_sd <= 0) {
    stop("initial_sd must be a positive finite number")
  }
  if (length(t) != length(y)) stop("t and y must be the same length")
  if (!is.numeric(t) || !is.numeric(y)) stop("t and y must be numeric vectors")

  t_obs_raw <- as.numeric(t)
  y_obs_raw <- as.numeric(y)

  # build a latent grid: unique observed times + infer_at
  time_grid_raw <- sort(unique(t_obs_raw))
  if (!is.null(infer_at)) {
    if (!is.numeric(infer_at)) stop("infer_at must be a numeric vector")
    time_grid_raw <- sort(unique(c(time_grid_raw, as.numeric(infer_at))))
  }

  if (length(time_grid_raw) < 2) {
    stop("Need at least 2 unique time points in t (or via infer_at) to define an IBM prior.")
  }

  # determine if grid is regular (using the grid bc determines trans density)
  regular <- all(abs(diff(time_grid_raw) - diff(time_grid_raw)[1]) == 0)

  # scale y to mean 0 sd 1
  y_mean <- mean(y_obs_raw)
  y_sd   <- stats::sd(y_obs_raw)
  y_obs  <- (y_obs_raw - y_mean) / y_sd

  # scale time using grid spacing so mean dt is 1ish
  deltat_raw <- diff(time_grid_raw)
  deltat_raw <- deltat_raw[deltat_raw > 0]
  dt_mean <- mean(deltat_raw)

  t_min <- min(time_grid_raw)
  time_grid <- (time_grid_raw - t_min) / dt_mean
  t_obs <- (t_obs_raw - t_min) / dt_mean

  # dt on SCALED grid not original
  deltat <- diff(time_grid)

  # map each observation to the grid index (replicates share the same index)
  obs_time_idx <- match(t_obs, time_grid)

  # make the stan data list
  stan_data <- list(
    N_obs = length(t_obs), # number of observations (can include replicates)
    T = length(time_grid), # number of latent grid points (unique + infer_at)
    obs_time_idx = obs_time_idx, # length N_obs, values in 1 to T
    deltat = deltat,  # length T-1, on scaled grid
    y_obs = y_obs, # length N_obs
    regular = as.integer(regular),
    initial_sd = initial_sd,
    log_sigma_mu = log_sigma$mu,
    log_sigma_sd = log_sigma$sd
  )

  if (adaptive %in% c("nonadaptive", "rw", "baseline_horseshoe", "persistent_clock")) {
    stan_data$log_tau_mu <- log_tau$mu
    stan_data$log_tau_sd <- log_tau$sd
  }
  if (adaptive %in% c("horseshoe", "baseline_horseshoe", "persistent_clock")) {
    stan_data$zeta <- zeta
  }
  if (adaptive == "rhs") {
    stan_data$zeta <- zeta
    stan_data$slab_scale <- slab_scale
  }
  if (adaptive == "bridge") {
    stan_data$alpha <- alpha
    stan_data$zeta <- zeta
  }
  # fit model with rstan
  stanfit <- rstan::stan(
    file = system.file(
      "stan",
      file,
      package = "ibmsmooth"
    ),
    data = stan_data,
    iter = iter,
    chains = chains,
    cores = cores,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta),
    ...
  )

  fit <- list(
    stanfit = stanfit,
    data = list(
      y_raw = y_obs_raw,
      y_mean = y_mean,
      y_sd = y_sd,
      t_raw = t_obs_raw,
      t_min = t_min,
      regular = regular,
      dt_mean = dt_mean,
      time_grid_raw = time_grid_raw,
      time_grid = time_grid,
      obs_time_idx = obs_time_idx,
      initial_sd = initial_sd
    ),
    adaptive = adaptive,
    centered = centered
  )
  class(fit) <- c("ibmfit", class(fit))
  return(fit)
}
