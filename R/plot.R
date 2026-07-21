#' Plot method for ibmfit objects
#'
#' Visualize fitted results from an object of class "ibmfit". The plot typically shows model
#' predictions with uncertainty and, if supplied, observed truth values for comparison.
#'
#' @param x An object of class "ibmfit", produced by fitting the IBM smoothing model.
#' @param truth Optional vector, data.frame, or time series of true/observed values to overlay
#'   on the plot. If provided, it will be used to compare fitted values to ground truth.
#' @param titles Optional character vector of length one or two giving plot titles for the function and derivative plots.
#' @param ... Additional graphical parameters passed to underlying plotting functions (e.g.,
#'   type, col, lwd) or other method-specific options recognized by plot.ibmfit.
#'
#' @details
#' This S3 method produces diagnostic and summary plots for an "ibmfit" object, such as fitted
#' trajectories, uncertainty intervals, and residual summaries. When `truth` is supplied,
#' observed values are added to the plot for visual comparison. The exact panels and layout
#' depend on the contents of the `ibmfit` object.
#'
#' @return A list with components:
#'   \item{function_plot}{A ggplot object showing the estimated function and its uncertainty.}
#'   \item{derivative_plot}{A ggplot object showing the estimated derivative and its uncertainty.}
#'
#' @import ggplot2
#'
#' @export
#' @method plot ibmfit
plot.ibmfit <- function(x, truth = NULL, titles = NULL, ...) {
  dat <- x$data

  # observed data (can have replicates)
  dat_orig <- data.frame(
    t = dat$t_raw,
    y = dat$y_raw
  )

  # latent grid in original time units (includes infer_at if used)
  # fallback for older fit objects
  if (!is.null(dat$time_grid_raw)) {
    t_grid_raw <- dat$time_grid_raw
  } else {
    t_grid_raw <- sort(unique(dat$t_raw))
  }

  # Determine whether this is a Stan or INLA fit
  if (!is.null(x$state_samples_scaled)) {
    f_samples_rescaled <- get_samples(x, "f")
    fprime_samples_rescaled <- get_samples(x, "fprime")
    if (is.null(titles)) {
      titles_use <- c("Marginalized horseshoe IBM: Estimated Function",
                      "Marginalized horseshoe IBM: Estimated Derivative")
    } else {
      titles_use <- titles
      if (length(titles_use) < 2L)
        titles_use <- c(titles_use[1L], "Marginalized horseshoe IBM: Estimated Derivative")
    }
    return(plot_curve(
      t_unique = t_grid_raw, f_samples = f_samples_rescaled,
      fprime_samples = fprime_samples_rescaled, dat_orig = dat_orig,
      truth = truth, titles = titles_use, ...
    ))
  }
  if (!is.null(x$stanfit)) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("The rstan package is required to plot a Stan fit.")
    }
    # Stan backend: extract posterior samples on latent grid
    stanfit <- x$stanfit
    f_samples <- rstan::extract(stanfit, pars = "f")$f
    fprime_samples <- rstan::extract(stanfit, pars = "fprime")$fprime
    # sanity check
    if (ncol(f_samples) != length(t_grid_raw)) {
      stop(
        "Mismatch: f has ", ncol(f_samples), " time points but grid has ",
        length(t_grid_raw), ". Did you change the Stan model/data conventions?"
      )
    }
    # rescale posterior samples back to original y units
    f_samples_rescaled <- f_samples * dat$y_sd + dat$y_mean
    # convert derivative to raw time scale
    fprime_samples_rescaled <- fprime_samples * (dat$y_sd / dat$dt_mean)
    # determine titles
    if (is.null(titles)) {
      titles_use <- c("IBM Smoothing: Estimated Function", "IBM Smoothing: Estimated Derivative")
    } else {
      titles_use <- titles
      if (length(titles_use) < 2) {
        titles_use <- c(titles_use[1], "IBM Smoothing: Estimated Derivative")
      }
    }
    return(plot_curve(
      t_unique = t_grid_raw,
      f_samples = f_samples_rescaled,
      fprime_samples = fprime_samples_rescaled,
      dat_orig = dat_orig,
      truth = truth,
      titles = titles_use,
      ...
    ))
  }
  if (!is.null(x$inla_obj)) {
    # INLA backend: sample from posterior using helper
    # by default draw a modest number of samples for plotting
    nsamp <- list(...)[["nsamp"]]
    if (is.null(nsamp)) nsamp <- 500
    f_samples <- get_samples_inla(x, param = "f", n_samples = nsamp)
    fprime_samples <- get_samples_inla(x, param = "fprime", n_samples = nsamp)
    # rescaling handled in get_samples_inla
    if (ncol(f_samples) != length(t_grid_raw)) {
      stop("Mismatch: sampled f has different length than latent grid")
    }
    if (is.null(titles)) {
      titles_use <- c("IBM Smoothing (INLA): Estimated Function", "IBM Smoothing (INLA): Estimated Derivative")
    } else {
      titles_use <- titles
      if (length(titles_use) < 2) {
        titles_use <- c(titles_use[1], "IBM Smoothing (INLA): Estimated Derivative")
      }
    }
    return(plot_curve(
      t_unique = t_grid_raw,
      f_samples = f_samples,
      fprime_samples = fprime_samples,
      dat_orig = dat_orig,
      truth = truth,
      titles = titles_use,
      ...
    ))
  }
  stop("Unrecognised ibmfit object: neither stanfit nor inla_obj found")
}

#' Plot posterior function and derivative summaries
#'
#' Construct ggplot objects showing posterior summaries for the latent function
#' and, optionally, its derivative.
#'
#' @param t_unique Numeric vector of time points on the latent grid.
#' @param f_samples Matrix of posterior samples for the latent function, with
#'   rows representing posterior draws and columns representing time points.
#' @param fprime_samples Optional matrix of posterior samples for the derivative,
#'   with rows representing posterior draws and columns representing time points.
#' @param dat_orig Optional data frame containing observed data with columns
#'   `t` and `y`.
#' @param truth Optional vector, list, or data frame containing truth values to
#'   overlay. Lists or data frames may contain `truth` and `deriv`.
#' @param titles Optional character vector of length one or two giving titles
#'   for the function and derivative plots.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing one or two ggplot objects.
#' @export
plot_curve <- function(t_unique, f_samples, fprime_samples = NULL,
                       dat_orig = NULL, truth = NULL, titles = NULL, ...) {

  # compute posterior summaries (95% and 80%)
  f_median  <- apply(f_samples, 2, stats::median)
  f_lower95 <- apply(f_samples, 2, stats::quantile, probs = 0.025)
  f_upper95 <- apply(f_samples, 2, stats::quantile, probs = 0.975)
  f_lower80 <- apply(f_samples, 2, stats::quantile, probs = 0.10)
  f_upper80 <- apply(f_samples, 2, stats::quantile, probs = 0.90)

  fit_df <- data.frame(
    t = t_unique,
    f_median = f_median,
    f_lower95 = f_lower95,
    f_upper95 = f_upper95,
    f_lower80 = f_lower80,
    f_upper80 = f_upper80
  )

  # derivative summaries only if provided
  if (!is.null(fprime_samples)) {
    fprime_median  <- apply(fprime_samples, 2, stats::median)
    fprime_lower95 <- apply(fprime_samples, 2, stats::quantile, probs = 0.025)
    fprime_upper95 <- apply(fprime_samples, 2, stats::quantile, probs = 0.975)
    fprime_lower80 <- apply(fprime_samples, 2, stats::quantile, probs = 0.10)
    fprime_upper80 <- apply(fprime_samples, 2, stats::quantile, probs = 0.90)

    fit_df$fprime_median  <- fprime_median
    fit_df$fprime_lower95 <- fprime_lower95
    fit_df$fprime_upper95 <- fprime_upper95
    fit_df$fprime_lower80 <- fprime_lower80
    fit_df$fprime_upper80 <- fprime_upper80
  }

  # flexible handling for truth input: must align to t_unique (i.e., latent grid)
  if (!is.null(truth)) {
    if (is.data.frame(truth)) {
      if (nrow(truth) != length(t_unique)) {
        stop("truth data.frame has ", nrow(truth), " rows but ", length(t_unique), " time points expected")
      }
      if ("truth" %in% names(truth)) fit_df$truth_f <- truth$truth
      if ("deriv" %in% names(truth)) fit_df$truth_fprime <- truth$deriv
    } else if (is.list(truth)) {
      if (!is.null(truth$truth) && length(truth$truth) != length(t_unique)) {
        stop("truth$truth has length ", length(truth$truth), " but ", length(t_unique), " time points expected")
      }
      fit_df$truth_f <- truth$truth
      if (!is.null(truth$deriv)) fit_df$truth_fprime <- truth$deriv
    } else if (is.numeric(truth)) {
      if (length(truth) != length(t_unique)) {
        stop("truth vector has length ", length(truth), " but ", length(t_unique), " time points expected")
      }
      fit_df$truth_f <- truth
    } else {
      stop("truth must be a data.frame, list, or numeric vector")
    }
  }

  # default titles if not provided
  if (is.null(titles)) {
    titles <- c("IBM Smoothing: Estimated Function", "IBM Smoothing: Estimated Derivative")
  } else if (length(titles) < 2) {
    titles <- c(titles[1], "IBM Smoothing: Estimated Derivative")
  }

  p1 <- ggplot2::ggplot() +
    {if (!is.null(dat_orig)) ggplot2::geom_point(data = dat_orig, ggplot2::aes(x = t, y = y), color = "gray") else NULL} +
    ggplot2::geom_ribbon(data = fit_df, ggplot2::aes(x = t, ymin = f_lower95, ymax = f_upper95), fill = "blue", alpha = 0.2) +
    ggplot2::geom_ribbon(data = fit_df, ggplot2::aes(x = t, ymin = f_lower80, ymax = f_upper80), fill = "blue", alpha = 0.2) +
    ggplot2::geom_line(data = fit_df, ggplot2::aes(x = t, y = f_median), color = "blue") +
    {if (!is.null(truth) && "truth_f" %in% names(fit_df))
      ggplot2::geom_line(data = fit_df, ggplot2::aes(x = t, y = truth_f), color = "black", linetype = "dashed")
    else NULL} +
    ggplot2::labs(title = titles[1], x = "t", y = "f(t)") +
    ggplot2::theme_minimal()

  plots <- list(function_plot = p1)

  if (!is.null(fprime_samples)) {
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = fit_df, ggplot2::aes(x = t, ymin = fprime_lower95, ymax = fprime_upper95), fill = "red", alpha = 0.2) +
      ggplot2::geom_ribbon(data = fit_df, ggplot2::aes(x = t, ymin = fprime_lower80, ymax = fprime_upper80), fill = "red", alpha = 0.2) +
      ggplot2::geom_line(data = fit_df, ggplot2::aes(x = t, y = fprime_median), color = "red") +
      {if (!is.null(truth) && "truth_fprime" %in% names(fit_df))
        ggplot2::geom_line(data = fit_df, ggplot2::aes(x = t, y = truth_fprime), color = "black", linetype = "dashed")
      else NULL} +
      ggplot2::labs(title = titles[2], x = "t", y = "f'(t)") +
      ggplot2::theme_minimal()

    plots$derivative_plot <- p2
  }

  plots
}

#' Plot the locally adaptive process precision curve
#'
#' Plot posterior summaries of the time-varying IBM process precision for a
#' locally adaptive fit. The curve is evaluated at transition midpoints, so each
#' plotted time value represents the interval between two adjacent latent grid
#' points.
#'
#' @param ibmfit An object of class `ibmfit` from a locally adaptive Stan or
#'   INLA fit.
#' @param n_samples Number of posterior samples to use. For INLA fits these are
#'   drawn with `INLA::inla.posterior.sample()`. For Stan fits this is the
#'   maximum number of MCMC draws to keep.
#' @param level Outer pointwise credible interval level. Defaults to `0.95`.
#' @param inner_level Optional inner pointwise credible interval level. Defaults
#'   to `0.80`. Use `NULL` to suppress the inner interval.
#' @param log_y Logical. If `TRUE`, display the positive precision values on a
#'   log10 y-axis.
#' @param title Optional plot title.
#' @param line_color Color for the posterior median curve.
#' @param ribbon_fill Fill color for the credible interval ribbons.
#' @param alpha_outer Alpha for the outer credible interval ribbon.
#' @param alpha_inner Alpha for the inner credible interval ribbon.
#' @param ... Additional arguments, currently unused.
#'
#' @return A `ggplot` object. The data used to construct the plot are attached
#'   as the `precision_summary` attribute.
#' @export
plot_precision_curve <- function(ibmfit, n_samples = 1000,
                                 level = 0.95,
                                 inner_level = 0.80,
                                 log_y = TRUE,
                                 title = NULL,
                                 line_color = "purple4",
                                 ribbon_fill = "purple",
                                 alpha_outer = 0.18,
                                 alpha_inner = 0.30,
                                 ...) {
  if (!inherits(ibmfit, "ibmfit")) {
    stop("ibmfit must be an object of class 'ibmfit'.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(inner_level) &&
      (!is.numeric(inner_level) || length(inner_level) != 1L ||
       !is.finite(inner_level) || inner_level <= 0 || inner_level >= 1)) {
    stop("inner_level must be NULL or a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(inner_level) && inner_level >= level) {
    stop("inner_level must be smaller than level.", call. = FALSE)
  }

  hp <- get_hyperparameter_samples(
    ibmfit,
    n_samples = n_samples,
    format = "long",
    natural = TRUE,
    include_internal = FALSE
  )

  if (nrow(hp) == 0L || !("process_precision" %in% hp$parameter)) {
    stop("No process-precision samples were found in this fit.", call. = FALSE)
  }

  precision_samples <- hp[hp$parameter == "process_precision" & !is.na(hp$t), , drop = FALSE]
  if (nrow(precision_samples) == 0L) {
    stop(
      "This fit appears to have a single global process precision. ",
      "plot_precision_curve() is intended for locally adaptive fits with ",
      "time-varying process precision.",
      call. = FALSE
    )
  }

  outer_lower_prob <- (1 - level) / 2
  outer_upper_prob <- 1 - outer_lower_prob
  if (is.null(inner_level)) {
    probs <- c(outer_lower_prob, 0.5, outer_upper_prob)
  } else {
    inner_lower_prob <- (1 - inner_level) / 2
    inner_upper_prob <- 1 - inner_lower_prob
    probs <- c(outer_lower_prob, inner_lower_prob, 0.5, inner_upper_prob, outer_upper_prob)
  }

  # Group by transition component. The index is stable for adaptive fits; fall
  # back to t in case older objects did not store an index.
  group_key <- ifelse(
    !is.na(precision_samples$index) & precision_samples$index != "",
    precision_samples$index,
    as.character(precision_samples$t)
  )
  groups <- split(precision_samples, group_key)

  plot_df <- do.call(rbind, lapply(groups, function(d) {
    q <- stats::quantile(d$value, probs = probs, na.rm = TRUE, names = FALSE)
    out <- data.frame(
      index = d$index[1L],
      t = d$t[1L],
      mean = mean(d$value, na.rm = TRUE),
      lower = q[1L],
      median = q[which.min(abs(probs - 0.5))],
      upper = q[length(q)],
      row.names = NULL
    )
    if (!is.null(inner_level)) {
      out$lower_inner <- q[which.min(abs(probs - inner_lower_prob))]
      out$upper_inner <- q[which.min(abs(probs - inner_upper_prob))]
    }
    out
  }))
  plot_df <- plot_df[order(plot_df$t), , drop = FALSE]

  if (is.null(title)) {
    title <- "Locally adaptive IBM process precision"
  }

  y_lab <- if (isTRUE(log_y)) "process precision (log10 scale)" else "process precision"

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = t, y = median)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = ribbon_fill,
      alpha = alpha_outer
    ) +
    {if (!is.null(inner_level))
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_inner, ymax = upper_inner),
        fill = ribbon_fill,
        alpha = alpha_inner
      )
    else NULL} +
    ggplot2::geom_line(color = line_color) +
    ggplot2::labs(
      title = title,
      x = "t",
      y = y_lab
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }

  attr(p, "precision_summary") <- plot_df
  p
}

#' Plot baseline-horseshoe IBM roughness components
#'
#' Plot posterior pointwise summaries of either total local roughness or the
#' horseshoe excess roughness. Values are process standard deviations on the
#' original response/time scale and are evaluated at transition midpoints.
#'
#' @param ibmfit A baseline-horseshoe `ibmfit` object.
#' @param component Either `"total"` or `"excess"`.
#' @param n_samples Maximum number of posterior draws to use.
#' @param level Pointwise credible interval level.
#' @param log_y Logical; use a log10 y-axis.
#' @param title Optional plot title.
#' @param line_color,ribbon_fill Plot colors.
#' @param ... Additional arguments, currently unused.
#'
#' @return A `ggplot` object with its summary data in the `roughness_summary`
#'   attribute.
#' @export
plot_roughness_curve <- function(ibmfit, component = c("total", "excess"),
                                 n_samples = 1000, level = 0.95,
                                 log_y = TRUE, title = NULL,
                                 line_color = "purple4",
                                 ribbon_fill = "purple", ...) {
  component <- match.arg(component)
  if (!inherits(ibmfit, "ibmfit")) {
    stop("ibmfit must be an object of class 'ibmfit'.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number between 0 and 1.", call. = FALSE)
  }
  parameter <- if (component == "total") "process_sd" else "excess_process_sd"
  hp <- get_hyperparameter_samples(ibmfit, n_samples = n_samples,
                                   format = "long", natural = TRUE)
  draws <- hp[hp$parameter == parameter & !is.na(hp$t), , drop = FALSE]
  if (nrow(draws) == 0L) {
    stop("No ", component, " local roughness samples were found in this fit.",
         call. = FALSE)
  }
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  groups <- split(draws, draws$index)
  plot_df <- do.call(rbind, lapply(groups, function(d) {
    q <- stats::quantile(d$value, probs = probs, na.rm = TRUE, names = FALSE)
    data.frame(index = d$index[1L], t = d$t[1L], mean = mean(d$value),
               lower = q[1L], median = q[2L], upper = q[3L])
  }))
  plot_df <- plot_df[order(plot_df$t), , drop = FALSE]
  if (is.null(title)) {
    title <- paste("Baseline-horseshoe IBM", component, "local roughness")
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = t, y = median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = ribbon_fill, alpha = 0.2) +
    ggplot2::geom_line(color = line_color) +
    ggplot2::labs(title = title, x = "t", y = paste(component, "process SD")) +
    ggplot2::theme_minimal()
  if (isTRUE(log_y)) p <- p + ggplot2::scale_y_log10()
  attr(p, "roughness_summary") <- plot_df
  p
}
