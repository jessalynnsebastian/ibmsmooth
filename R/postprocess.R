#' Summarize an ibmfit object
#'
#' Compute pointwise posterior summaries for the latent function, its
#' derivative, and the available model hyperparameters. Unlike a scalar
#' parameter summary, the curve summaries are returned at every point of the
#' latent time grid.
#'
#' @param object An object of class `"ibmfit"`.
#' @param ... Additional arguments. Currently passed to curve-summary methods.
#' @param n_samples Number of posterior samples to draw for INLA summaries. For
#'   Stan fits, existing MCMC draws are used; if `n_samples` is smaller than the
#'   number of available draws, the first `n_samples` draws are used.
#' @param probs Numeric vector of posterior probabilities to summarize.
#'
#' @return A list with components `curve` and `hyperparameters`.
#' @export
#' @method summary ibmfit
summary.ibmfit <- function(object, ..., n_samples = 1000,
                           probs = c(0.025, 0.10, 0.50, 0.90, 0.975)) {
  out <- list(
    curve = get_curve_summary(object, n_samples = n_samples, probs = probs),
    hyperparameters = get_hyperparameter_summary(object, n_samples = n_samples, probs = probs)
  )
  class(out) <- "summary_ibmfit"
  out
}

#' Print an ibmfit summary
#'
#' @param x An object returned by `summary.ibmfit()`.
#' @param ... Additional arguments, currently unused.
#'
#' @return The input object, invisibly.
#' @export
#' @method print summary_ibmfit
print.summary_ibmfit <- function(x, ...) {
  cat("ibmfit posterior summary\n")
  cat("\nCurve summaries are pointwise over the latent grid. Use `$curve` to access all rows.\n")
  print(utils::head(x$curve, 10))
  if (!is.null(x$hyperparameters) && nrow(x$hyperparameters) > 0) {
    cat("\nHyperparameter summaries. Use `$hyperparameters` to access all rows.\n")
    print(x$hyperparameters)
  } else {
    cat("\nNo hyperparameter summaries available.\n")
  }
  invisible(x)
}

#' Extract posterior curve samples
#'
#' Extract posterior samples for the latent function and/or derivative over the
#' full latent grid.
#'
#' @param ibmfit An object of class `ibmfit`.
#' @param param Character. One of `"f"`, `"fprime"`, or `"both"`.
#' @param n_samples Number of samples to use for INLA fits, or maximum number
#'   of MCMC draws to keep for Stan fits. Use `NULL` to keep all Stan draws.
#' @param format Output format, either `"matrix"` or `"long"`.
#'
#' @return If `format = "matrix"` and one parameter is requested, a posterior
#'   sample matrix with rows as draws and columns as time points. If both
#'   parameters are requested, a list of two matrices. If `format = "long"`, a
#'   data frame with columns `draw`, `t`, `parameter`, and `value`.
#' @export
get_curve_samples <- function(ibmfit, param = c("f", "fprime", "both"),
                              n_samples = 1000,
                              format = c("matrix", "long")) {
  param <- match.arg(param)
  format <- match.arg(format)
  dat <- ibmfit$data
  t_grid <- dat$time_grid_raw

  params <- if (param == "both") c("f", "fprime") else param
  mats <- stats::setNames(vector("list", length(params)), params)

  for (p in params) {
    mats[[p]] <- get_samples(ibmfit, param = p, n_samples = n_samples)
  }

  if (format == "matrix") {
    if (length(mats) == 1L) return(mats[[1L]])
    return(mats)
  }

  pieces <- lapply(names(mats), function(p) {
    mat <- mats[[p]]
    data.frame(
      draw = rep(seq_len(nrow(mat)), each = ncol(mat)),
      t = rep(t_grid, times = nrow(mat)),
      parameter = p,
      value = as.numeric(t(mat)),
      row.names = NULL
    )
  })
  do.call(rbind, pieces)
}

#' Summarize posterior curves
#'
#' Compute pointwise posterior summaries for the latent function and derivative.
#'
#' @param ibmfit An object of class `ibmfit`.
#' @param n_samples Number of samples to use for INLA fits, or maximum number
#'   of MCMC draws to keep for Stan fits.
#' @param probs Numeric vector of posterior probabilities.
#'
#' @return A data frame with one row per parameter and time point.
#' @export
get_curve_summary <- function(ibmfit, n_samples = 1000,
                              probs = c(0.025, 0.10, 0.50, 0.90, 0.975)) {
  dat <- ibmfit$data
  t_grid <- dat$time_grid_raw
  mats <- get_curve_samples(ibmfit, param = "both", n_samples = n_samples, format = "matrix")

  pieces <- lapply(names(mats), function(p) {
    mat <- mats[[p]]
    q <- apply(mat, 2, stats::quantile, probs = probs, na.rm = TRUE)
    if (is.vector(q)) q <- matrix(q, nrow = length(probs))
    q_df <- as.data.frame(t(q), check.names = FALSE)
    names(q_df) <- paste0("q", formatC(100 * probs, format = "f", digits = 1))
    data.frame(
      t = t_grid,
      parameter = p,
      mean = colMeans(mat, na.rm = TRUE),
      sd = apply(mat, 2, stats::sd, na.rm = TRUE),
      q_df,
      row.names = NULL,
      check.names = FALSE
    )
  })
  do.call(rbind, pieces)
}

#' Extract posterior hyperparameter samples
#'
#' Extract posterior samples for model hyperparameters. By default, returned
#' quantities are transformed to the original response and time scales.
#'
#' For INLA fits, hyperparameters are sampled with `intern = TRUE`, so the
#' Gaussian likelihood hyperparameter is log precision and the rgeneric field
#' hyperparameters are log diffusion precisions on the scaled model scale. On
#' the natural scale this function
#' returns the implied IBM process standard deviation and precision. For
#' adaptive INLA fits these are returned over transition midpoints, and the
#' global baseline process standard deviation and precision are also returned
#' as `global_process_sd` and `global_process_precision`.
#'
#' For Stan fits, recognized transformed parameters such as `sigma`, `tau`, and
#' `gamma` are converted from scaled model units. Baseline-horseshoe fits also
#' report `baseline_process_sd`, `global_excess_sd`, `excess_process_sd`,
#' `process_sd` (total roughness), and `process_variance`.
#'
#' @param ibmfit An object of class `ibmfit`.
#' @param n_samples Number of samples to draw for INLA fits, or maximum number
#'   of Stan draws to keep. Use `NULL` to keep all Stan draws.
#' @param format Output format, either `"long"` or `"wide"`.
#' @param natural Logical. If `TRUE`, return quantities transformed to the
#'   original response/time scale whenever possible. If `FALSE`, return internal
#'   backend hyperparameters.
#' @param include_internal Logical. If `TRUE` and `natural = TRUE`, also include
#'   internal backend hyperparameters.
#'
#' @return A data frame of hyperparameter samples.
#' @export
get_hyperparameter_samples <- function(ibmfit, n_samples = 1000,
                                       format = c("long", "wide"),
                                       natural = TRUE,
                                       include_internal = FALSE) {
  format <- match.arg(format)
  dat <- ibmfit$data

  make_long_from_array <- function(arr, name, keep, value_transform = identity,
                                   parameter = name, t = NA_real_) {
    if (is.null(dim(arr)) || length(dim(arr)) == 1L) {
      arr <- matrix(arr, ncol = 1)
    }
    arr <- arr[keep, , drop = FALSE]
    trailing_dim <- dim(arr)[-1]
    mat <- matrix(arr, nrow = length(keep))
    if (ncol(mat) == 1L) {
      index <- ""
      t_vec <- t[1L]
    } else {
      index <- apply(expand.grid(lapply(trailing_dim, seq_len)), 1, paste, collapse = ",")
      t_vec <- if (length(t) == ncol(mat)) t else rep(NA_real_, ncol(mat))
    }
    mat <- value_transform(mat)
    do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
      data.frame(
        draw = keep,
        parameter = parameter,
        index = index[j],
        t = t_vec[j],
        value = mat[, j],
        row.names = NULL
      )
    }))
  }

  if (!is.null(ibmfit$stanfit)) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("The rstan package is required to extract hyperparameter samples from a Stan fit.")
    }
    samples <- rstan::extract(ibmfit$stanfit, permuted = TRUE)
    samples <- samples[setdiff(names(samples), c("f", "fprime", "lp__"))]
    if (length(samples) == 0L) return(data.frame())

    n_draw <- dim(samples[[1]])[1]
    if (is.null(n_draw)) n_draw <- length(samples[[1]])
    keep <- seq_len(n_draw)
    if (!is.null(n_samples) && n_draw > n_samples) keep <- seq_len(n_samples)

    pieces <- list()
    process_scale <- dat$y_sd / (dat$dt_mean^(3 / 2))

    if (isTRUE(natural)) {
      if (!is.null(samples$sigma)) {
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$sigma, "sigma", keep,
          value_transform = function(x) x * dat$y_sd,
          parameter = "observation_sd"
        )
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$sigma, "sigma", keep,
          value_transform = function(x) 1 / ((x * dat$y_sd)^2),
          parameter = "observation_precision"
        )
      }
      if (!is.null(samples$tau)) {
        transition_t <- ibmfit$data$time_grid_raw[-1]
        if (!is.null(ibmfit$inla_model$transition_midpoints_raw)) {
          transition_t <- ibmfit$inla_model$transition_midpoints_raw
        } else if (length(dat$time_grid_raw) > 1L) {
          transition_t <- (dat$time_grid_raw[-1] + dat$time_grid_raw[-length(dat$time_grid_raw)]) / 2
        }
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$tau, "tau", keep,
          value_transform = function(x) x * process_scale,
          parameter = "process_sd",
          t = transition_t
        )
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$tau, "tau", keep,
          value_transform = function(x) 1 / ((x * process_scale)^2),
          parameter = "process_precision",
          t = transition_t
        )
      }
      if (!is.null(samples$tau0)) {
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$tau0, "tau0", keep,
          value_transform = function(x) x * process_scale,
          parameter = "baseline_process_sd"
        )
      }
      if (!is.null(samples$psi)) {
        transition_t <- (dat$time_grid_raw[-1] +
          dat$time_grid_raw[-length(dat$time_grid_raw)]) / 2
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$psi, "psi", keep,
          value_transform = function(x) x * process_scale,
          parameter = "excess_process_sd", t = transition_t
        )
      }
      if (!is.null(samples$local_variance)) {
        transition_t <- (dat$time_grid_raw[-1] +
          dat$time_grid_raw[-length(dat$time_grid_raw)]) / 2
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$local_variance, "local_variance", keep,
          value_transform = function(x) x * process_scale^2,
          parameter = "process_variance", t = transition_t
        )
      }
      if (!is.null(samples$gamma)) {
        pieces[[length(pieces) + 1L]] <- make_long_from_array(
          samples$gamma, "gamma", keep,
          value_transform = function(x) x * process_scale,
          parameter = if (!is.null(samples$tau0)) "global_excess_sd" else "global_process_sd"
        )
      }
    }

    if (!isTRUE(natural) || isTRUE(include_internal)) {
      internal <- do.call(rbind, lapply(names(samples), function(nm) {
        make_long_from_array(samples[[nm]], nm, keep, parameter = paste0("internal_", nm))
      }))
      pieces[[length(pieces) + 1L]] <- internal
    }

    long <- if (length(pieces) == 0L) data.frame() else do.call(rbind, pieces)

  } else if (!is.null(ibmfit$inla_obj)) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("The INLA package is required to sample from an INLA fit.")
    }
    # Important: by default, INLA::inla.posterior.sample() returns likelihood
    # hyperparameters on their external/natural INLA scale. For the Gaussian
    # likelihood that means precision, not log precision. The transformations
    # below are written on the internal theta scale, so request intern = TRUE
    # explicitly. This keeps the Gaussian likelihood hyperparameter as log
    # precision and leaves rgeneric theta parameters on the scale used by Q().
    intern <- TRUE
    samples_list <- tryCatch(
      INLA::inla.posterior.sample(n = n_samples, result = ibmfit$inla_obj, intern = TRUE),
      error = function(e) {
        if (grepl("unused argument", conditionMessage(e), fixed = TRUE)) {
          intern <<- FALSE
          INLA::inla.posterior.sample(n = n_samples, result = ibmfit$inla_obj)
        } else {
          stop(
            "INLA posterior sampling failed. Refit the model with ",
            "`control.compute = list(config = TRUE)`. Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      }
    )
    if (length(samples_list) == 0L || is.null(samples_list[[1]]$hyperpar)) return(data.frame())

    hp <- do.call(rbind, lapply(samples_list, function(s) as.numeric(s$hyperpar)))
    colnames(hp) <- names(samples_list[[1]]$hyperpar)
    keep <- seq_len(nrow(hp))

    pieces <- list()

    obs_col <- grep("Log precision for the Gaussian observations", colnames(hp), fixed = TRUE)
    if (length(obs_col) == 0L) {
      obs_col <- grep("Gaussian observations", colnames(hp), fixed = TRUE)
    }
    field_cols <- grep("for field", colnames(hp), fixed = TRUE)

    if (isTRUE(natural)) {
      if (length(obs_col) >= 1L) {
        if (isTRUE(intern)) {
          obs_log_prec <- hp[, obs_col[1L]]
          obs_prec_scaled <- exp(obs_log_prec)
        } else {
          # Older INLA fallback: default posterior samples return Gaussian
          # precision on INLA's external scale.
          obs_prec_scaled <- hp[, obs_col[1L]]
        }
        obs_sd <- dat$y_sd / sqrt(obs_prec_scaled)
        obs_prec <- obs_prec_scaled / (dat$y_sd^2)
        pieces[[length(pieces) + 1L]] <- data.frame(
          draw = keep, parameter = "observation_sd", index = "",
          t = NA_real_, value = obs_sd, row.names = NULL
        )
        pieces[[length(pieces) + 1L]] <- data.frame(
          draw = keep, parameter = "observation_precision", index = "",
          t = NA_real_, value = obs_prec, row.names = NULL
        )
      }

      if (length(field_cols) > 0L) {
        theta <- hp[, field_cols, drop = FALSE]
        B <- ibmfit$inla_model$B
        if (is.null(B)) {
          stop("This INLA fit does not store the B matrix needed to transform field hyperparameters. Please refit the model.")
        }
        if (ncol(theta) != ncol(B)) {
          stop("Could not match INLA field hyperparameters to the stored B matrix.")
        }

        log_prec_scaled <- theta %*% t(B)
        process_sd <- exp(-0.5 * log_prec_scaled) * dat$y_sd / (dat$dt_mean^(3 / 2))
        process_precision <- 1 / (process_sd^2)
        transition_t <- ibmfit$inla_model$transition_midpoints_raw
        if (is.null(transition_t)) {
          transition_t <- (dat$time_grid_raw[-1] + dat$time_grid_raw[-length(dat$time_grid_raw)]) / 2
        }

        # If this is the nonadaptive INLA model, report the global process
        # hyperparameter once rather than repeating it at every interval.
        is_global <- ncol(B) == 1L && length(unique(as.numeric(B))) == 1L
        is_global_plus_deviation <- identical(
          ibmfit$inla_model$adaptive_parameterization,
          "global_plus_spline_deviation"
        )

        if (is_global_plus_deviation) {
          global_log_prec_scaled <- theta[, 1L]
          global_process_sd <- exp(-0.5 * global_log_prec_scaled) *
            dat$y_sd / (dat$dt_mean^(3 / 2))
          global_process_precision <- 1 / (global_process_sd^2)

          pieces[[length(pieces) + 1L]] <- data.frame(
            draw = keep, parameter = "global_process_sd", index = "",
            t = NA_real_, value = as.numeric(global_process_sd), row.names = NULL
          )
          pieces[[length(pieces) + 1L]] <- data.frame(
            draw = keep, parameter = "global_process_precision", index = "",
            t = NA_real_, value = as.numeric(global_process_precision), row.names = NULL
          )
        }

        if (is_global) {
          pieces[[length(pieces) + 1L]] <- data.frame(
            draw = keep, parameter = "process_sd", index = "",
            t = NA_real_, value = as.numeric(process_sd[, 1L]), row.names = NULL
          )
          pieces[[length(pieces) + 1L]] <- data.frame(
            draw = keep, parameter = "process_precision", index = "",
            t = NA_real_, value = as.numeric(process_precision[, 1L]), row.names = NULL
          )
        } else {
          pieces[[length(pieces) + 1L]] <- do.call(rbind, lapply(seq_len(ncol(process_sd)), function(j) {
            data.frame(
              draw = keep, parameter = "process_sd", index = as.character(j),
              t = transition_t[j], value = process_sd[, j], row.names = NULL
            )
          }))
          pieces[[length(pieces) + 1L]] <- do.call(rbind, lapply(seq_len(ncol(process_precision)), function(j) {
            data.frame(
              draw = keep, parameter = "process_precision", index = as.character(j),
              t = transition_t[j], value = process_precision[, j], row.names = NULL
            )
          }))
        }
      }
    }

    if (!isTRUE(natural) || isTRUE(include_internal)) {
      internal_prefix <- if (isTRUE(intern)) "internal_" else "external_inla_"
      internal <- do.call(rbind, lapply(seq_len(ncol(hp)), function(j) {
        data.frame(
          draw = keep,
          parameter = paste0(internal_prefix, colnames(hp)[j]),
          index = "",
          t = NA_real_,
          value = hp[, j],
          row.names = NULL
        )
      }))
      pieces[[length(pieces) + 1L]] <- internal
    }

    long <- if (length(pieces) == 0L) data.frame() else do.call(rbind, pieces)

  } else {
    stop("Unrecognised ibmfit object: cannot extract hyperparameter samples")
  }

  if (format == "long") return(long)
  if (nrow(long) == 0L) return(data.frame())

  long$key <- ifelse(
    is.na(long$t),
    ifelse(long$index == "", long$parameter, paste0(long$parameter, "[", long$index, "]")),
    paste0(long$parameter, "[", long$index, "]")
  )
  stats::reshape(
    long[, c("draw", "key", "value")],
    idvar = "draw",
    timevar = "key",
    direction = "wide"
  )
}

#' Summarize posterior hyperparameters
#'
#' @param ibmfit An object of class `ibmfit`.
#' @param n_samples Number of samples to draw for INLA fits, or maximum number
#'   of Stan draws to keep.
#' @param probs Numeric vector of posterior probabilities.
#' @param natural Logical. If `TRUE`, summarize hyperparameters on the original
#'   response/time scale when possible.
#' @param include_internal Logical. If `TRUE`, include internal backend
#'   hyperparameters in addition to natural-scale quantities.
#'
#' @return A data frame with one row per hyperparameter component.
#' @export
get_hyperparameter_summary <- function(ibmfit, n_samples = 1000,
                                       probs = c(0.025, 0.10, 0.50, 0.90, 0.975),
                                       natural = TRUE,
                                       include_internal = FALSE) {
  hp <- get_hyperparameter_samples(
    ibmfit,
    n_samples = n_samples,
    format = "long",
    natural = natural,
    include_internal = include_internal
  )
  if (nrow(hp) == 0L) return(data.frame())

  combos <- unique(hp[, c("parameter", "index", "t"), drop = FALSE])
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    parameter_i <- combos$parameter[i]
    index_i <- combos$index[i]
    t_i <- combos$t[i]
    idx_match <- hp$parameter == parameter_i & hp$index == index_i
    if (is.na(t_i)) {
      idx_match <- idx_match & is.na(hp$t)
    } else {
      idx_match <- idx_match & !is.na(hp$t) & hp$t == t_i
    }
    values <- hp$value[idx_match]
    q <- stats::quantile(values, probs = probs, na.rm = TRUE)
    q_df <- as.data.frame(as.list(q), check.names = FALSE)
    names(q_df) <- paste0("q", formatC(100 * probs, format = "f", digits = 1))
    data.frame(
      parameter = parameter_i,
      index = index_i,
      t = t_i,
      mean = mean(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      q_df,
      row.names = NULL,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Extract posterior samples from an IBM fit
#'
#' @param ibmfit An object of class `ibmfit`.
#' @param param Character, either `"f"` or `"fprime"`.
#' @param n_samples Optional number of posterior draws to return.
#' @param ... Additional arguments passed to `get_samples_inla()` for INLA fits.
#'
#' @return A matrix of posterior samples.
#' @export
get_samples <- function(ibmfit, param = c("f", "B_operational", "fprime_left",
                                          "fprime_right", "fprime"), n_samples = NULL, ...) {
  param <- match.arg(param)

  if (!is.null(ibmfit$state_samples_scaled)) {
    samples <- ibmfit$state_samples_scaled[[param]]
    dat <- ibmfit$data
    samples_rescaled <- if (param == "f") {
      samples * dat$y_sd + dat$y_mean
    } else {
      samples * (dat$y_sd / dat$dt_mean)
    }
    if (!is.null(n_samples) && nrow(samples_rescaled) > n_samples)
      samples_rescaled <- samples_rescaled[seq_len(n_samples), , drop = FALSE]
    return(samples_rescaled)
  }

  if (!is.null(ibmfit$stanfit)) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("The rstan package is required to extract samples from a Stan fit.")
    }
    stan_param <- param
    if (param == "fprime" && ibmfit$adaptive %in% c("nonadaptive", "baseline_horseshoe", "persistent_clock")) {
      warning("`fprime` is deprecated for stochastic-clock fits; returning right derivatives with the final knot undefined.")
      stan_param <- "fprime_right"
    }
    samples <- rstan::extract(ibmfit$stanfit, pars = stan_param)[[stan_param]]
    dat <- ibmfit$data
    if (param == "f") {
      samples_rescaled <- samples * dat$y_sd + dat$y_mean
    } else if (param %in% c("fprime", "fprime_left", "fprime_right")) {
      samples_rescaled <- samples * (dat$y_sd / dat$dt_mean)
    } else {
      samples_rescaled <- samples * dat$y_sd
    }
    if (param == "fprime_left") samples_rescaled[, 1L] <- NA_real_
    if (param %in% c("fprime", "fprime_right") &&
        ibmfit$adaptive %in% c("nonadaptive", "baseline_horseshoe", "persistent_clock")) {
      samples_rescaled[, ncol(samples_rescaled)] <- NA_real_
    }
    if (!is.null(n_samples) && nrow(samples_rescaled) > n_samples) {
      samples_rescaled <- samples_rescaled[seq_len(n_samples), , drop = FALSE]
    }
    return(samples_rescaled)
  }

  if (!is.null(ibmfit$inla_obj)) {
    return(get_samples_inla(
      ibmfit,
      param = param,
      n_samples = if (is.null(n_samples)) 1000 else n_samples,
      ...
    ))
  }

  stop("Unrecognised ibmfit object: cannot extract samples")
}

.ibmsmooth_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

.ibmsmooth_inla_latent_index <- function(result, latent_names, field_name, n_state) {
  escaped <- .ibmsmooth_escape_regex(field_name)

  patterns <- c(
    paste0("^", escaped, "\\["),
    paste0("^", escaped, ":"),
    paste0("^", escaped, "\\."),
    paste0("^", escaped, "\\("),
    paste0("^", escaped, "$")
  )

  idx <- unique(unlist(lapply(patterns, grep, x = latent_names)))
  if (length(idx) == n_state) return(idx)

  contents <- tryCatch(result$misc$configs$contents, error = function(e) NULL)

  if (is.data.frame(contents) && all(c("tag", "start", "length") %in% names(contents))) {
    row <- which(contents$tag == field_name)
    if (length(row) == 1L && contents$length[row] == n_state) {
      return(seq.int(contents$start[row] + 1L, length.out = n_state))
    }
    rows <- which(contents$tag != "Predictor" & contents$length == n_state)
    if (length(rows) == 1L) {
      return(seq.int(contents$start[rows] + 1L, length.out = n_state))
    }
  }

  if (is.list(contents)) {
    if (!is.null(contents$tag) && !is.null(contents$start) && !is.null(contents$length)) {
      row <- which(contents$tag == field_name)
      if (length(row) == 1L && contents$length[row] == n_state) {
        return(seq.int(contents$start[row] + 1L, length.out = n_state))
      }
    }
    if (!is.null(contents[[field_name]])) {
      entry <- contents[[field_name]]
      if (!is.null(entry$start) && !is.null(entry$length) && entry$length == n_state) {
        return(seq.int(entry$start + 1L, length.out = n_state))
      }
    }
  }

  stop(
    "Could not identify the INLA rgeneric latent field in posterior samples. ",
    "First latent names were: ",
    paste(utils::head(latent_names, 10), collapse = ", "),
    call. = FALSE
  )
}

#' Extract posterior samples from an INLA fit
#'
#' Helper function to draw posterior samples for the latent function or its
#' derivative from an INLA fit. Samples are drawn using
#' `INLA::inla.posterior.sample()`. The latent field is rescaled back to
#' the original observation scale.
#'
#' @param ibmfit An object of class `ibmfit` produced by the INLA implementation.
#' @param param Character, either `"f"` for the latent function or `"fprime"` for its derivative.
#' @param n_samples Integer giving the number of posterior samples to draw.
#'
#' @return A matrix of dimension `n_samples` by length of the latent time grid.
#' @export
get_samples_inla <- function(ibmfit, param = c("f", "fprime"), n_samples = 1000) {
  param <- match.arg(param, c("f", "fprime"))
  if (is.null(ibmfit$inla_obj)) stop("Input does not contain an INLA fit")
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required to sample from an INLA fit.")
  }

  res <- ibmfit$inla_obj
  dat <- ibmfit$data
  n_grid <- length(dat$time_grid)
  n_state <- 2L * n_grid

  samples_list <- tryCatch(
    INLA::inla.posterior.sample(n = n_samples, result = res),
    error = function(e) {
      stop(
        "INLA posterior sampling failed. Refit the model with ",
        "`control.compute = list(config = TRUE)`. This is now the default ",
        "for `ibm_inla_fit()` and `ibm(method = \"inla\")`. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  latent_names <- rownames(samples_list[[1]]$latent)

  idx_latent <- .ibmsmooth_inla_latent_index(
    result = res,
    latent_names = latent_names,
    field_name = ibmfit$rgeneric_name,
    n_state = n_state
  )

  out <- matrix(NA_real_, nrow = n_samples, ncol = n_grid)
  for (i in seq_len(n_samples)) {
    latent <- as.numeric(samples_list[[i]]$latent[idx_latent])
    if (param == "f") {
      out[i, ] <- latent[seq_len(n_grid)] * dat$y_sd + dat$y_mean
    } else {
      out[i, ] <- latent[n_grid + seq_len(n_grid)] * (dat$y_sd / dat$dt_mean)
    }
  }
  out
}
