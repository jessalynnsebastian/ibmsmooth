#' Fit an integrated Brownian motion smoother using INLA
#'
#' Fit a global or locally adaptive integrated Brownian motion smoother using
#' INLA. The latent field is the state vector `(f, fprime)` on the latent time
#' grid. The observation model is Gaussian.
#'
#' @param t Numeric vector of time points.
#' @param y Numeric vector of observations at times `t`.
#' @param infer_at Optional numeric vector of additional time points at which
#'   inference should be returned.
#' @param adaptive Logical. If `FALSE`, use one global diffusion-precision parameter. If
#'   `TRUE`, model the log diffusion precision as a global log precision plus
#'   centered B-spline deviations shrunk toward zero.
#' @param n_knots Number of internal knots for the B-spline basis when
#'   `adaptive = TRUE`.
#' @param bs_degree Degree of the B-spline basis.
#' @param adaptive_deviation_sd Prior standard deviation for the centered
#'   B-spline deviation coefficients in adaptive INLA fits. The adaptive model
#'   uses `log(lambda_i) = alpha + B_i beta`, where `alpha` is the global log
#'   diffusion precision and the components of `beta` have independent
#'   `N(0, adaptive_deviation_sd^2)` priors. Smaller values shrink the model
#'   more strongly toward the global-precision IBM.
#' @param log_process_sd List with elements `mu` and `sd` giving the prior
#'   mean and standard deviation for the log IBM process standard deviation on
#'   the scaled model scale. The default matches the Stan backend:
#'   `list(mu = -2, sd = 0.5)`.
#' @param log_sigma List with elements `mu` and `sd` giving the prior mean and
#'   standard deviation for the log Gaussian observation standard deviation on
#'   the scaled model scale. The default matches the Stan backend:
#'   `list(mu = -1, sd = 1)`.
#' @param theta_mu Optional advanced override for the mean of the internal INLA
#'   global log diffusion-precision parameter. If `NULL`, this is set from
#'   `log_process_sd`.
#' @param theta_sd Optional advanced override for the standard deviation of the
#'   internal INLA global log diffusion-precision parameter. If `NULL`, this is
#'   set from `log_process_sd`.
#' @param sigma_mu Optional advanced override for the mean of the Gaussian
#'   likelihood log precision. If `NULL`, this is set from `log_sigma`.
#' @param sigma_sd Optional advanced override for the standard deviation of the
#'   Gaussian likelihood log precision. If `NULL`, this is set from `log_sigma`.
#' @param initial_sd Prior scale for the joint initial state `c(f_1, fprime_1)`.
#'   The initial covariance is `initial_sd^2 * matrix(c(1/3, 1/2, 1/2, 1), 2, 2)`,
#'   the unit-time integrated Brownian motion covariance for the level and derivative.
#' @param initial_slope_sd Deprecated. Retained only for backward compatibility;
#'   ignored by the INLA backend.
#' @param control.compute List passed to `INLA::inla()` as `control.compute`.
#'   Defaults to `list(config = TRUE)` so that posterior sampling works for
#'   plotting, summaries, and `get_samples_inla()`.
#' @param ... Additional arguments passed to `INLA::inla()`.
#'
#' @return An object of class `ibmfit`.
#' @export
ibm_inla_fit <- function(t, y, infer_at = NULL, adaptive = FALSE,
                         n_knots = 5L, bs_degree = 3L,
                         adaptive_deviation_sd = 0.5,
                         log_process_sd = list(mu = -2, sd = 0.5),
                         log_sigma = list(mu = -1, sd = 1),
                         theta_mu = NULL, theta_sd = NULL,
                         sigma_mu = NULL, sigma_sd = NULL,
                         initial_sd = 5, initial_slope_sd = NULL,
                         control.compute = list(config = TRUE), ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required for method = 'inla'. Install it from the INLA repository.")
  }
  if (length(t) != length(y)) stop("t and y must be the same length")
  if (!is.numeric(t) || !is.numeric(y)) stop("t and y must be numeric vectors")
  if (!is.numeric(initial_sd) || length(initial_sd) != 1 ||
      !is.finite(initial_sd) || initial_sd <= 0) {
    stop("initial_sd must be a positive finite number")
  }
  # `initial_slope_sd` is retained only so older code does not break. The INLA
  # backend uses one correlated initial-state covariance controlled by
  # `initial_sd`, so this argument is ignored.
  if (!is.null(initial_slope_sd)) {
    initial_slope_sd <- NULL
  }

  if (is.null(log_process_sd$mu) || is.null(log_process_sd$sd) || log_process_sd$sd <= 0) {
    stop("log_process_sd must be a list with numeric elements mu and positive sd")
  }
  if (is.null(log_sigma$mu) || is.null(log_sigma$sd) || log_sigma$sd <= 0) {
    stop("log_sigma must be a list with numeric elements mu and positive sd")
  }

  # INLA uses log precisions internally. Convert natural log-standard-deviation
  # priors to log-precision priors:
  #   log_precision = -2 * log_sd.
  # Thus log_sd ~ N(mu, sd^2) implies log_precision ~ N(-2mu, (2sd)^2).
  if (is.null(theta_mu)) theta_mu <- -2 * log_process_sd$mu
  if (is.null(theta_sd)) theta_sd <- 2 * log_process_sd$sd
  if (is.null(sigma_mu)) sigma_mu <- -2 * log_sigma$mu
  if (is.null(sigma_sd)) sigma_sd <- 2 * log_sigma$sd

  if (!is.numeric(theta_mu) || length(theta_mu) != 1 || !is.finite(theta_mu)) {
    stop("theta_mu must be a finite numeric scalar")
  }
  if (!is.numeric(theta_sd) || length(theta_sd) != 1 || !is.finite(theta_sd) || theta_sd <= 0) {
    stop("theta_sd must be a positive finite numeric scalar")
  }
  if (!is.numeric(adaptive_deviation_sd) || length(adaptive_deviation_sd) != 1 ||
      !is.finite(adaptive_deviation_sd) || adaptive_deviation_sd <= 0) {
    stop("adaptive_deviation_sd must be a positive finite numeric scalar")
  }
  if (!is.numeric(sigma_mu) || length(sigma_mu) != 1 || !is.finite(sigma_mu)) {
    stop("sigma_mu must be a finite numeric scalar")
  }
  if (!is.numeric(sigma_sd) || length(sigma_sd) != 1 || !is.finite(sigma_sd) || sigma_sd <= 0) {
    stop("sigma_sd must be a positive finite numeric scalar")
  }

  if (is.null(control.compute)) {
    control.compute <- list(config = TRUE)
  }
  if (is.null(control.compute$config)) {
    control.compute$config <- TRUE
  }

  t_obs_raw <- as.numeric(t)
  y_obs_raw <- as.numeric(y)

  time_grid_raw <- sort(unique(t_obs_raw))
  if (!is.null(infer_at)) {
    if (!is.numeric(infer_at)) stop("infer_at must be numeric")
    time_grid_raw <- sort(unique(c(time_grid_raw, as.numeric(infer_at))))
  }
  if (length(time_grid_raw) < 2) {
    stop("Need at least 2 unique time points in t or infer_at.")
  }

  regular <- all(abs(diff(time_grid_raw) - diff(time_grid_raw)[1]) < sqrt(.Machine$double.eps))

  y_mean <- mean(y_obs_raw)
  y_sd <- stats::sd(y_obs_raw)
  if (!is.finite(y_sd) || y_sd <= 0) stop("y must have positive finite standard deviation")
  y_obs <- (y_obs_raw - y_mean) / y_sd

  dt_raw <- diff(time_grid_raw)
  dt_mean <- mean(dt_raw[dt_raw > 0])
  t_min <- min(time_grid_raw)
  time_grid <- (time_grid_raw - t_min) / dt_mean
  t_obs <- (t_obs_raw - t_min) / dt_mean
  dt <- diff(time_grid)

  n_grid <- length(time_grid)
  n_state <- 2L * n_grid
  obs_time_idx <- match(t_obs, time_grid)

  transition_midpoints <- (time_grid[-1] + time_grid[-n_grid]) / 2
  transition_midpoints_raw <- (time_grid_raw[-1] + time_grid_raw[-n_grid]) / 2

  if (isTRUE(adaptive)) {
    n_knots <- as.integer(n_knots)
    if (n_knots < 0) stop("n_knots must be non-negative")
    internal_knots <- NULL
    if (n_knots > 0) {
      probs <- seq(0, 1, length.out = n_knots + 2L)
      internal_knots <- as.numeric(stats::quantile(transition_midpoints, probs = probs))
      internal_knots <- internal_knots[seq.int(2L, length(internal_knots) - 1L)]
    }
    B_raw <- splines::bs(
      transition_midpoints,
      knots = internal_knots,
      degree = bs_degree,
      intercept = TRUE
    )
    B_raw <- as.matrix(B_raw)

    # Adaptive INLA parameterization:
    #   log(lambda_i) = alpha + B_dev[i, ] %*% beta,
    # where alpha is the global log diffusion precision and beta is shrunk
    # toward zero.  Center the spline columns so beta represents departures
    # from the global precision rather than absorbing the overall level.
    B_centered <- sweep(B_raw, 2L, colMeans(B_raw), FUN = "-")
    if (ncol(B_centered) > 1L) {
      # After centering, the full intercept B-spline basis has one redundant
      # column because the original rows sum to one. Drop one column so INLA
      # does not integrate over a likelihood-null deviation direction.
      B_deviation <- B_centered[, seq_len(ncol(B_centered) - 1L), drop = FALSE]
    } else {
      B_deviation <- matrix(0, nrow = n_grid - 1L, ncol = 0L)
    }
    B <- cbind(global = rep(1, n_grid - 1L), B_deviation)
    if (ncol(B) > 1L) {
      colnames(B) <- c("global", paste0("deviation", seq_len(ncol(B) - 1L)))
    }
  } else {
    B_raw <- NULL
    B_deviation <- NULL
    B <- matrix(1, nrow = n_grid - 1L, ncol = 1L)
    colnames(B) <- "global"
  }
  nbasis <- ncol(B)
  theta_mu_vec <- c(theta_mu, rep(0, nbasis - 1L))
  theta_sd_vec <- c(theta_sd, rep(adaptive_deviation_sd, nbasis - 1L))

  # Joint correlated prior for the initial state c(f_1, fprime_1), on the
  # scaled response and scaled time scales. The base covariance is the
  # unit-time IBM covariance for the level and derivative.
  initial_cov <- initial_sd^2 * matrix(c(1/3, 1/2, 1/2, 1), nrow = 2)
  initial_prec <- solve(initial_cov)
  log_det_initial_prec <- as.numeric(determinant(initial_prec, logarithm = TRUE)$modulus)
  log_det_transition_base <- log(12) - 4 * log(dt)

  # INLA rgeneric models are evaluated in an INLA-controlled context.
  # To avoid scoping failures, embed all fit-specific objects directly into
  # the rgeneric function body using substitute().  This makes the callback
  # self-contained when INLA serializes/evaluates it internally.
  ibm_rgeneric <- eval(substitute(
    function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
             theta = NULL) {
      cmd <- match.arg(cmd)

      B <- B_VALUE
      dt <- DT_VALUE
      n_grid <- N_GRID_VALUE
      n_state <- N_STATE_VALUE
      nbasis <- NBASIS_VALUE
      initial_prec <- INITIAL_PREC_VALUE
      log_det_initial_prec <- LOG_DET_INITIAL_PREC_VALUE
      log_det_transition_base <- LOG_DET_TRANSITION_BASE_VALUE
      theta_mu_vec <- THETA_MU_VALUE
      theta_sd_vec <- THETA_SD_VALUE

      if (!length(theta)) theta <- theta_mu_vec

      make_precision_local <- function(theta) {
        eta <- as.numeric(B %*% theta)
        trans_prec <- exp(eta)

        # Joint correlated initial-state prior precision for
        # c(f_1, fprime_1). With symmetric = TRUE, Matrix expects only one
        # triangle, so use the upper-triangular entries of the 2 x 2 block.
        initial_idx <- c(1L, n_grid + 1L)
        Q <- Matrix::sparseMatrix(
          i = c(initial_idx[1L], initial_idx[1L], initial_idx[2L]),
          j = c(initial_idx[1L], initial_idx[2L], initial_idx[2L]),
          x = c(initial_prec[1L, 1L], initial_prec[1L, 2L], initial_prec[2L, 2L]),
          dims = c(n_state, n_state),
          symmetric = TRUE
        )

        # Transition precision for the IBM state:
        #   f'_i - f'_{i-1}
        #   f_i  - f_{i-1} - h f'_{i-1}
        # with covariance proportional to
        #   [[h, h^2 / 2], [h^2 / 2, h^3 / 3]].
        for (i in 2:n_grid) {
          h <- dt[i - 1L]
          # Error vector is ordered as:
          #   e_1 = fprime_i - fprime_{i-1}
          #   e_2 = f_i - f_{i-1} - h * fprime_{i-1}
          # so the IBM innovation covariance is for c(e_1, e_2):
          #   [[h, h^2 / 2], [h^2 / 2, h^3 / 3]].
          S <- matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), nrow = 2)
          Sinv <- solve(S)

          H <- Matrix::sparseMatrix(
            i = c(1L, 1L, 2L, 2L, 2L),
            j = c(n_grid + i, n_grid + i - 1L, i, i - 1L, n_grid + i - 1L),
            x = c(1, -1, 1, -1, -h),
            dims = c(2L, n_state)
          )

          transition_q <- trans_prec[i - 1L] *
            Matrix::t(H) %*% Matrix::Matrix(Sinv, sparse = TRUE) %*% H

          # Keep Q explicitly symmetric. This avoids Matrix validity issues when
          # sparse products contain both triangles or tiny numerical asymmetries.
          Q <- Q + Matrix::forceSymmetric(transition_q, uplo = "U")
        }

        Matrix::drop0(Matrix::forceSymmetric(Q, uplo = "U"))
      }

      log_norm_const_local <- function(theta) {
        if (!length(theta)) theta <- theta_mu_vec
        eta <- as.numeric(B %*% theta)

        # The rgeneric density contribution is -0.5 * x'Qx plus this
        # normalizing constant. Compute it analytically from the product of
        # the initial-state density and the IBM transition densities. This is
        # less error-prone than taking a determinant of the assembled sparse
        # precision and makes the diffusion-precision scaling explicit.
        #
        # For one transition, with theta_i = log(lambda_i),
        #   Cov(e_i | theta_i) = lambda_i^{-1}
        #     [[h, h^2 / 2], [h^2 / 2, h^3 / 3]],
        # where e_i = (fprime_i - fprime_{i-1},
        #              f_i - f_{i-1} - h fprime_{i-1}).
        # Therefore 0.5 * log |Q_i| = theta_i +
        # 0.5 * log |S_i^{-1}|, and |S_i^{-1}| = 12 / h^4.
        0.5 * log_det_initial_prec +
          sum(eta + 0.5 * log_det_transition_base) -
          0.5 * n_state * log(2 * pi)
      }

      if (cmd == "graph") {
        return(make_precision_local(theta_mu_vec))
      }
      if (cmd == "Q") {
        return(make_precision_local(theta))
      }
      if (cmd == "mu") {
        return(numeric(n_state))
      }
      if (cmd == "initial") {
        return(theta_mu_vec)
      }
      if (cmd == "log.norm.const") {
        return(log_norm_const_local(theta))
      }
      if (cmd == "log.prior") {
        return(sum(stats::dnorm(theta, mean = theta_mu_vec, sd = theta_sd_vec, log = TRUE)))
      }
      if (cmd == "quit") {
        return(invisible())
      }

      stop("Unknown rgeneric command")
    },
    list(
      B_VALUE = B,
      DT_VALUE = dt,
      N_GRID_VALUE = n_grid,
      N_STATE_VALUE = n_state,
      NBASIS_VALUE = nbasis,
      INITIAL_PREC_VALUE = initial_prec,
      LOG_DET_INITIAL_PREC_VALUE = log_det_initial_prec,
      LOG_DET_TRANSITION_BASE_VALUE = log_det_transition_base,
      THETA_MU_VALUE = theta_mu_vec,
      THETA_SD_VALUE = theta_sd_vec
    )
  ))

  # Do not leave the function environment pointing at ibm_inla_fit's frame:
  # everything needed by the callback has been embedded in the function body.
  environment(ibm_rgeneric) <- baseenv()

  rgeneric_model <- INLA::inla.rgeneric.define(
    model = ibm_rgeneric,
    n = n_state
  )

  # Observations use only the level component of the state vector.  The full
  # latent state has length 2T, ordered as f[1:T], fprime[1:T].
  dat_inla <- data.frame(y = y_obs, field = obs_time_idx)
  formula_inla <- y ~ -1 + f(field, model = rgeneric_model, values = seq_len(n_state))

  # INLA's Gaussian likelihood hyperparameter is log precision. A normal prior
  # uses parameters c(mean, precision), so the second entry is 1 / sd^2.
  hyper_prec <- list(prec = list(prior = "normal", param = c(sigma_mu, 1 / sigma_sd^2)))

  result <- INLA::inla(
    formula_inla,
    data = dat_inla,
    family = "gaussian",
    control.family = list(hyper = hyper_prec),
    control.predictor = list(compute = TRUE),
    control.compute = control.compute,
    ...
  )

  fit <- list(
    inla_obj = result,
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
      n_state = n_state
    ),
    rgeneric_name = "field",
    inla_model = list(
      adaptive = adaptive,
      n_knots = n_knots,
      bs_degree = bs_degree,
      initial_sd = initial_sd,
      initial_cov = initial_cov,
      B = B,
      B_raw = B_raw,
      B_deviation = B_deviation,
      nbasis = nbasis,
      adaptive_parameterization = if (isTRUE(adaptive)) "global_plus_spline_deviation" else "global",
      transition_midpoints = transition_midpoints,
      transition_midpoints_raw = transition_midpoints_raw,
      transition_parameter = "diffusion_precision",
      transition_error_order = c("fprime_increment", "f_increment"),
      log_process_sd = log_process_sd,
      log_sigma = log_sigma,
      theta_mu = theta_mu,
      theta_sd = theta_sd,
      theta_mu_vec = theta_mu_vec,
      theta_sd_vec = theta_sd_vec,
      adaptive_deviation_sd = adaptive_deviation_sd,
      sigma_mu = sigma_mu,
      sigma_sd = sigma_sd
    )
  )
  class(fit) <- c("ibmfit", class(fit))
  fit
}
