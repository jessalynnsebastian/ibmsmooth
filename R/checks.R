#' Draw from the prior associated with a fitted model
#'
#' The prior fit uses the same locations, scaling, and prior settings as
#' `ibmfit`, but omits the observation likelihood.
#'
#' @param ibmfit An `ibmfit` object returned by [ibm()].
#' @param iter,chains,cores Stan sampling controls for the prior fit.
#' @param ... Additional arguments passed to `rstan::sampling()`, such as
#'   `seed` or `refresh`.
#' @return An `ibmfit` containing prior draws.
#' @export
sample_prior <- function(ibmfit, iter = 1000, chains = 2,
                         cores = getOption("mc.cores", chains), ...) {
  .check_ibmfit(ibmfit)
  spec <- ibmfit$fit_spec
  if (is.null(spec)) {
    stop("This fit predates stored model specifications; refit it first.",
         call. = FALSE)
  }
  do.call(
    .ibm_fit,
    c(spec, list(
      iter = iter, chains = chains, cores = cores, get_code = FALSE,
      prior_only = TRUE
    ), list(...))
  )
}

#' Plot prior and posterior distributions
#'
#' @param ibmfit A posterior `ibmfit`.
#' @param prior_fit Optional prior-only fit from [sample_prior()]. When omitted,
#'   a prior fit is sampled automatically.
#' @param parameters Hyperparameters to display. The default shows `sigma` and
#'   the model's global smoothing scale (`tau` or `gamma`).
#' @param n_samples Maximum draws from each distribution.
#' @param ... Arguments passed to [sample_prior()] when `prior_fit` is omitted.
#' @return A faceted ggplot object.
#' @export
plot_prior_posterior <- function(ibmfit, prior_fit = NULL, parameters = NULL,
                                 n_samples = 1000, ...) {
  .check_ibmfit(ibmfit)
  if (is.null(prior_fit)) prior_fit <- sample_prior(ibmfit, ...)
  .check_ibmfit(prior_fit)
  if (!isTRUE(prior_fit$prior_only)) {
    warning("prior_fit is not marked as a prior-only fit.", call. = FALSE)
  }
  if (is.null(parameters)) {
    parameters <- c("sigma", if (isTRUE(ibmfit$adaptive)) "gamma" else "tau")
  }
  if (!is.character(parameters) || !length(parameters)) {
    stop("parameters must be a non-empty character vector.", call. = FALSE)
  }
  collect <- function(fit, distribution) {
    do.call(rbind, lapply(parameters, function(parameter) {
      draws <- get_samples(fit, parameter, n_samples)
      data.frame(
        parameter = if (ncol(draws) == 1L) parameter else
          paste0(parameter, "[", rep(seq_len(ncol(draws)),
                                     each = nrow(draws)), "]"),
        distribution = distribution, value = as.numeric(draws),
        row.names = NULL
      )
    }))
  }
  d <- rbind(collect(prior_fit, "Prior"), collect(ibmfit, "Posterior"))
  d$distribution <- factor(d$distribution, c("Prior", "Posterior"))
  ggplot2::ggplot(
    d, ggplot2::aes(x = value, colour = distribution, fill = distribution)
  ) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.8,
                          na.rm = TRUE) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::scale_colour_manual(values = c(Prior = "grey40",
                                             Posterior = "#0072B2")) +
    ggplot2::scale_fill_manual(values = c(Prior = "grey70",
                                           Posterior = "#56B4E9")) +
    ggplot2::labs(x = NULL, y = "Density", colour = NULL, fill = NULL,
                  title = "Prior and posterior distributions") +
    ggplot2::theme_minimal()
}

#' Generate posterior predictive observations
#'
#' @param ibmfit A posterior `ibmfit`.
#' @param n_samples Number of posterior predictive data sets.
#' @param seed Optional random seed.
#' @return A matrix with predictive data sets in rows and observations in
#'   columns, in the original response units.
#' @export
posterior_predict <- function(ibmfit, n_samples = 500, seed = NULL) {
  .check_ibmfit(ibmfit)
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive number.", call. = FALSE)
  }
  n_samples <- as.integer(n_samples)
  f <- get_samples(ibmfit, "f", n_samples)
  sigma <- get_samples(ibmfit, "sigma", n_samples)[, 1L]
  n_samples <- nrow(f)
  mu <- f[, ibmfit$data$obs_time_idx, drop = FALSE]
  if (!is.null(seed)) set.seed(seed)
  mu + matrix(
    stats::rnorm(length(mu), sd = rep(sigma, ncol(mu))),
    nrow = n_samples
  )
}

#' Visual posterior predictive check
#'
#' @param ibmfit A posterior `ibmfit`.
#' @param type Either `"interval"` for observation-level predictive intervals
#'   or `"density"` for observed and replicated-data densities.
#' @param n_samples Number of replicated data sets.
#' @param probs Predictive interval probabilities for `type = "interval"`.
#' @param n_density Number of replicated densities to draw.
#' @param seed Optional random seed.
#' @return A ggplot object.
#' @export
plot_pp_check <- function(ibmfit, type = c("interval", "density"),
                          n_samples = 500, probs = c(0.05, 0.5, 0.95),
                          n_density = 30, seed = NULL) {
  type <- match.arg(type)
  yrep <- posterior_predict(ibmfit, n_samples, seed)
  observed <- ibmfit$data$y_raw
  if (type == "interval") {
    if (length(probs) != 3L || any(!is.finite(probs)) ||
        any(probs < 0 | probs > 1) || any(diff(probs) <= 0)) {
      stop("probs must contain three increasing probabilities in [0, 1].",
           call. = FALSE)
    }
    qs <- apply(yrep, 2L, stats::quantile, probs = probs)
    d <- data.frame(
      t = ibmfit$data$t_raw, observed = observed,
      lower = qs[1L, ], center = qs[2L, ], upper = qs[3L, ]
    )
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = t)) +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = lower, ymax = upper),
          colour = "#56B4E9", alpha = 0.7
        ) +
        ggplot2::geom_point(ggplot2::aes(y = center), colour = "#0072B2",
                            shape = 1) +
        ggplot2::geom_point(ggplot2::aes(y = observed), colour = "grey20",
                            size = 1.8) +
        ggplot2::labs(title = "Posterior predictive check",
                      subtitle = "Observed points and predictive intervals",
                      x = "t", y = "y") +
        ggplot2::theme_minimal()
    )
  }
  n_density <- min(as.integer(n_density), nrow(yrep))
  replicated <- data.frame(
    value = as.numeric(t(yrep[seq_len(n_density), , drop = FALSE])),
    draw = factor(rep(seq_len(n_density), each = ncol(yrep)))
  )
  ggplot2::ggplot(replicated, ggplot2::aes(x = value, group = draw)) +
    ggplot2::geom_density(colour = "#56B4E9", alpha = 0.35,
                          linewidth = 0.45) +
    ggplot2::geom_density(
      data = data.frame(value = observed),
      ggplot2::aes(x = value), inherit.aes = FALSE,
      colour = "grey15", linewidth = 1.1
    ) +
    ggplot2::labs(title = "Posterior predictive check",
                  subtitle = "Observed density (dark) and replicated densities",
                  x = "y", y = "Density") +
    ggplot2::theme_minimal()
}

.check_ibmfit <- function(x) {
  if (!inherits(x, "ibmfit") || is.null(x$stanfit) || is.null(x$data)) {
    stop("ibmfit must contain a Stan fit.", call. = FALSE)
  }
  invisible(x)
}
