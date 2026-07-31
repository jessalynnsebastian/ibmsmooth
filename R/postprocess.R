.stan_draws <- function(ibmfit, param, n_samples = NULL) {
  if (!inherits(ibmfit, "ibmfit") || is.null(ibmfit$stanfit)) {
    stop("ibmfit must contain a Stan fit.", call. = FALSE)
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required to extract draws.", call. = FALSE)
  }
  out <- rstan::extract(ibmfit$stanfit, pars = param,
                        permuted = TRUE)[[param]]
  if (is.null(out)) stop("Parameter not found: ", param, call. = FALSE)
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  if (!is.null(n_samples) && nrow(out) > n_samples) {
    out <- out[seq_len(n_samples), , drop = FALSE]
  }
  out
}

#' Extract posterior samples
#'
#' @param ibmfit An `ibmfit` object.
#' @param param Stan parameter name. Common choices are `"f"`, `"fprime"`,
#'   `"sigma"`, `"tau"`, `"gamma"`, `"xi"`, `"xi_regularized"`, `"slab"`,
#'   and `"lambda_interval"`.
#' @param n_samples Maximum number of draws; `NULL` keeps all draws.
#' @return A matrix with posterior draws in rows.
#' @export
get_samples <- function(ibmfit, param = c("f", "fprime"),
                        n_samples = NULL) {
  param <- match.arg(
    param,
    c("f", "fprime", "sigma", "tau", "gamma", "xi",
      "xi_regularized", "slab", "tau_interval", "lambda_interval")
  )
  draws <- .stan_draws(ibmfit, param, n_samples)
  dat <- ibmfit$data
  if (param == "f") {
    draws <- draws * dat$y_sd + dat$y_mean
  } else if (param == "fprime") {
    draws <- draws * dat$y_sd / dat$dt_mean
  } else if (param == "sigma") {
    draws <- draws * dat$y_sd
  } else if (param %in% c("tau", "gamma", "slab", "tau_interval")) {
    draws <- draws * dat$y_sd / dat$dt_mean^(3 / 2)
  } else if (param == "lambda_interval") {
    draws <- draws * dat$y_sd^2 / dat$dt_mean^3
  }
  draws
}

#' Extract posterior curve samples
#'
#' @param ibmfit An `ibmfit` object.
#' @param param One of `"f"`, `"fprime"`, or `"both"`.
#' @param n_samples Maximum number of draws.
#' @param format Either `"matrix"` or `"long"`.
#' @param new_t Optional finite locations at which to extract either curve.
#' @param seed Optional seed for conditional bridge draws at new locations.
#' @return A matrix, list of matrices, or long data frame.
#' @export
get_curve_samples <- function(ibmfit, param = c("f", "fprime", "both"),
                              n_samples = 1000,
                              format = c("matrix", "long"),
                              new_t = NULL, seed = NULL) {
  param <- match.arg(param)
  format <- match.arg(format)
  requested <- if (param == "both") c("f", "fprime") else param
  prediction <- predict_curve(
    ibmfit, new_t = new_t, n_samples = n_samples, seed = seed
  )
  draws <- prediction[requested]
  if (format == "matrix") {
    if (length(draws) == 1L) return(draws[[1L]])
    return(draws)
  }
  do.call(rbind, lapply(names(draws), function(name) {
    x <- draws[[name]]
    data.frame(
      draw = rep(seq_len(nrow(x)), each = ncol(x)),
      t = rep(prediction$t, nrow(x)),
      parameter = name,
      value = as.numeric(t(x)),
      row.names = NULL
    )
  }))
}

#' Summarize posterior curves
#'
#' @param ibmfit An `ibmfit` object.
#' @param n_samples Maximum number of draws.
#' @param probs Posterior probabilities.
#' @param param One of `"f"`, `"fprime"`, or `"both"`.
#' @param new_t Optional finite locations at which to summarize either curve.
#' @param seed Optional seed for conditional bridge draws at new locations.
#' @return A data frame with pointwise summaries.
#' @export
get_curve_summary <- function(ibmfit, n_samples = 1000,
                              probs = c(0.025, 0.5, 0.975),
                              param = c("both", "f", "fprime"),
                              new_t = NULL, seed = NULL) {
  param <- match.arg(param)
  draws <- get_curve_samples(
    ibmfit, param, n_samples, "matrix", new_t = new_t, seed = seed
  )
  if (param != "both") draws <- setNames(list(draws), param)
  prediction_t <- as.numeric(colnames(draws[[1L]]))
  do.call(rbind, lapply(names(draws), function(name) {
    x <- draws[[name]]
    q <- t(apply(x, 2L, stats::quantile, probs = probs))
    colnames(q) <- paste0("q", formatC(100 * probs, format = "f", digits = 1))
    data.frame(
      t = prediction_t, parameter = name,
      mean = colMeans(x), sd = apply(x, 2L, stats::sd),
      q, row.names = NULL, check.names = FALSE
    )
  }))
}

.evaluate_curve_truth <- function(truth, t, name) {
  if (is.function(truth)) truth <- truth(t)
  if (!is.numeric(truth) || length(truth) != length(t) ||
      any(!is.finite(truth))) {
    stop(
      "truth$", name, " must be a function of t or a finite numeric vector ",
      "with one value per evaluation location.", call. = FALSE
    )
  }
  as.numeric(truth)
}

.empirical_crps <- function(draws, truth) {
  n <- nrow(draws)
  first_term <- colMeans(abs(sweep(draws, 2L, truth, "-")))
  if (n == 1L) return(first_term)
  coefficients <- 2 * seq_len(n) - n - 1
  pair_term <- vapply(seq_len(ncol(draws)), function(j) {
    sum(coefficients * sort(draws[, j])) / n^2
  }, numeric(1L))
  first_term - pair_term
}

.curve_performance <- function(draws, truth, t, level, parameter) {
  alpha <- 1 - level
  lower <- apply(draws, 2L, stats::quantile, probs = alpha / 2,
                 names = FALSE)
  upper <- apply(draws, 2L, stats::quantile, probs = 1 - alpha / 2,
                 names = FALSE)
  estimate <- apply(draws, 2L, stats::median)
  pointwise <- data.frame(
    t = t, parameter = parameter, truth = truth, estimate = estimate,
    absolute_deviation = abs(estimate - truth),
    interval_width = upper - lower,
    covered = truth >= lower & truth <= upper,
    crps = .empirical_crps(draws, truth),
    lower = lower, upper = upper, row.names = NULL
  )
  aggregate <- data.frame(
    parameter = parameter,
    MAD = mean(pointwise$absolute_deviation),
    MCIW = mean(pointwise$interval_width),
    Coverage = 100 * mean(pointwise$covered),
    CRPS = mean(pointwise$crps),
    level = level, n_locations = length(t), row.names = NULL
  )
  list(aggregate = aggregate, pointwise = pointwise)
}

#' Summarize function and derivative estimation performance
#'
#' Computes mean absolute deviation (MAD) of the posterior median, mean
#' credible interval width (MCIW), empirical coverage percentage, and the
#' empirical continuous ranked probability score (CRPS), pointwise and averaged
#' across evaluation locations.
#'
#' @param ibmfit An `ibmfit` object.
#' @param truth Named list containing `f` and/or `fprime`. Each element can be a
#'   numeric vector evaluated at `new_t` or a function that accepts `new_t`.
#' @param new_t Optional finite evaluation locations. By default, uses the
#'   fitted prediction grid.
#' @param level Credible interval confidence level strictly between zero and
#'   one.
#' @param n_samples Maximum number of posterior draws.
#' @param seed Optional seed for conditional bridge draws at new locations.
#' @param pointwise If `TRUE`, return both aggregate and location-specific
#'   results; otherwise return only the aggregate data frame.
#' @return A data frame with one row per supplied truth, or a list containing
#'   `aggregate` and `pointwise` data frames when `pointwise = TRUE`.
#' @export
summarize_performance <- function(ibmfit, truth, new_t = NULL, level = 0.95,
                                  n_samples = 1000, seed = NULL,
                                  pointwise = FALSE) {
  if (!is.list(truth) || is.null(names(truth)) ||
      !any(c("f", "fprime") %in% names(truth))) {
    stop("truth must be a named list containing f and/or fprime.",
         call. = FALSE)
  }
  unknown <- setdiff(names(truth), c("f", "fprime"))
  if (length(unknown)) {
    stop("Unknown truth component: ", unknown[1L], ".", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a finite scalar strictly between 0 and 1.",
         call. = FALSE)
  }
  prediction <- predict_curve(
    ibmfit, new_t = new_t, n_samples = n_samples, seed = seed
  )
  requested <- intersect(c("f", "fprime"), names(truth))
  results <- lapply(requested, function(name) {
    truth_values <- .evaluate_curve_truth(truth[[name]], prediction$t, name)
    .curve_performance(
      prediction[[name]], truth_values, prediction$t, level, name
    )
  })
  aggregate <- do.call(rbind, lapply(results, `[[`, "aggregate"))
  rownames(aggregate) <- NULL
  if (!isTRUE(pointwise)) return(aggregate)
  pointwise_result <- do.call(rbind, lapply(results, `[[`, "pointwise"))
  rownames(pointwise_result) <- NULL
  list(aggregate = aggregate, pointwise = pointwise_result)
}

#' Extract posterior hyperparameter samples
#'
#' @param ibmfit An `ibmfit` object.
#' @param n_samples Maximum number of draws.
#' @param format Either `"long"` or `"wide"`.
#' @return A data frame.
#' @export
get_hyperparameter_samples <- function(ibmfit, n_samples = 1000,
                                       format = c("long", "wide")) {
  format <- match.arg(format)
  names <- c("sigma", if (isTRUE(ibmfit$adaptive)) "gamma" else "tau")
  if (isTRUE(ibmfit$adaptive) &&
      identical(ibmfit$adaptive_prior, "regularized_horseshoe")) {
    names <- c(names, "slab")
  }
  if (isTRUE(ibmfit$adaptive)) names <- c(names, "lambda_interval")
  pieces <- lapply(names, function(name) {
    x <- get_samples(ibmfit, name, n_samples)
    if (ncol(x) == 1L) {
      return(data.frame(draw = seq_len(nrow(x)), parameter = name,
                        index = NA_integer_, t = NA_real_, value = x[, 1L]))
    }
    midpoint <- (ibmfit$data$time_grid_raw[-1L] +
                   ibmfit$data$time_grid_raw[-length(ibmfit$data$time_grid_raw)]) / 2
    do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
      data.frame(draw = seq_len(nrow(x)), parameter = name,
                 index = j, t = midpoint[j], value = x[, j])
    }))
  })
  long <- do.call(rbind, pieces)
  if (format == "long") return(long)
  long$key <- ifelse(is.na(long$index), long$parameter,
                     paste0(long$parameter, "[", long$index, "]"))
  stats::reshape(long[c("draw", "key", "value")], idvar = "draw",
                 timevar = "key", direction = "wide")
}

#' Summarize posterior hyperparameters
#'
#' @inheritParams get_hyperparameter_samples
#' @param probs Posterior probabilities.
#' @return A data frame.
#' @export
get_hyperparameter_summary <- function(ibmfit, n_samples = 1000,
                                       probs = c(0.025, 0.5, 0.975)) {
  x <- get_hyperparameter_samples(ibmfit, n_samples, "long")
  group_index <- ifelse(is.na(x$index), "scalar", as.character(x$index))
  groups <- split(x, interaction(x$parameter, group_index, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    q <- stats::quantile(d$value, probs = probs)
    data.frame(
      parameter = d$parameter[1L], index = d$index[1L], t = d$t[1L],
      mean = mean(d$value), sd = stats::sd(d$value),
      as.list(q), row.names = NULL, check.names = FALSE
    )
  }))
}

#' Summarize an IBM fit
#'
#' @param object An `ibmfit` object.
#' @param ... Unused.
#' @param n_samples Maximum number of draws.
#' @param probs Posterior probabilities.
#' @return An object of class `summary_ibmfit`.
#' @export
#' @method summary ibmfit
summary.ibmfit <- function(object, ..., n_samples = 1000,
                           probs = c(0.025, 0.5, 0.975)) {
  out <- list(
    curve = get_curve_summary(object, n_samples, probs),
    hyperparameters = get_hyperparameter_summary(object, n_samples, probs)
  )
  class(out) <- "summary_ibmfit"
  out
}

#' Print an IBM fit summary
#'
#' @param x A `summary_ibmfit` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @method print summary_ibmfit
print.summary_ibmfit <- function(x, ...) {
  cat("ibmfit posterior summary\n\n")
  print(utils::head(x$curve, 10L))
  cat("\nHyperparameters\n")
  print(x$hyperparameters)
  invisible(x)
}
