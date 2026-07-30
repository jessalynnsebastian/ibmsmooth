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
#' @return A matrix, list of matrices, or long data frame.
#' @export
get_curve_samples <- function(ibmfit, param = c("f", "fprime", "both"),
                              n_samples = 1000,
                              format = c("matrix", "long")) {
  param <- match.arg(param)
  format <- match.arg(format)
  requested <- if (param == "both") c("f", "fprime") else param
  prediction <- predict_curve(ibmfit, n_samples = n_samples)
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
#' @return A data frame with pointwise summaries.
#' @export
get_curve_summary <- function(ibmfit, n_samples = 1000,
                              probs = c(0.025, 0.5, 0.975)) {
  draws <- get_curve_samples(ibmfit, "both", n_samples)
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
