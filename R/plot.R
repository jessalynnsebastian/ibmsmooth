#' Plot an integrated Brownian motion fit
#'
#' @param x An `ibmfit` object.
#' @param truth Optional list with `f` and/or `fprime` values on the latent grid.
#' @param ... Additional arguments passed to [plot_curve()].
#' @return A list containing function and derivative ggplots.
#' @export
#' @method plot ibmfit
plot.ibmfit <- function(x, truth = NULL, ...) {
  if (!inherits(x, "ibmfit")) stop("x must be an ibmfit object.", call. = FALSE)
  prediction_grid <- x$data$prediction_grid_raw
  if (is.null(prediction_grid)) prediction_grid <- x$data$time_grid_raw
  curves <- predict_curve(
    x, new_t = prediction_grid
  )
  plot_curve(
    t_unique = curves$t,
    f_samples = curves$f,
    fprime_samples = curves$fprime,
    dat_orig = data.frame(t = x$data$t_raw, y = x$data$y_raw),
    truth = truth,
    ...
  )
}

#' Plot posterior function and derivative summaries
#'
#' @param t_unique Latent-grid locations.
#' @param f_samples Posterior function draws, with draws in rows.
#' @param fprime_samples Optional posterior derivative draws.
#' @param dat_orig Optional data frame with `t` and `y`.
#' @param truth Optional list with `f` and/or `fprime`.
#' @param titles Optional two plot titles.
#' @param ... Unused.
#' @return A list of ggplot objects.
#' @export
plot_curve <- function(t_unique, f_samples, fprime_samples = NULL,
                       dat_orig = NULL, truth = NULL, titles = NULL, ...) {
  summarize <- function(x) {
    data.frame(
      median = apply(x, 2L, stats::median),
      lower = apply(x, 2L, stats::quantile, probs = 0.025),
      upper = apply(x, 2L, stats::quantile, probs = 0.975)
    )
  }
  fsum <- cbind(t = t_unique, summarize(f_samples))
  if (is.null(titles)) {
    titles <- c("Estimated function", "Estimated derivative")
  }
  p1 <- ggplot2::ggplot(fsum, ggplot2::aes(x = t, y = median)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#56B4E9", alpha = 0.25
    ) +
    ggplot2::geom_line(colour = "#0072B2") +
    {if (!is.null(dat_orig))
      ggplot2::geom_point(
        data = dat_orig, ggplot2::aes(x = t, y = y),
        inherit.aes = FALSE, colour = "grey40"
      )
    else NULL} +
    {if (!is.null(truth) && !is.null(truth$f))
      ggplot2::geom_line(
        data = data.frame(t = t_unique, value = truth$f),
        ggplot2::aes(x = t, y = value), inherit.aes = FALSE,
        linetype = 2
      )
    else NULL} +
    ggplot2::labs(title = titles[1L], x = "t", y = "f(t)") +
    ggplot2::theme_minimal()
  out <- list(function_plot = p1)

  if (!is.null(fprime_samples)) {
    dsum <- cbind(t = t_unique, summarize(fprime_samples))
    out$derivative_plot <-
      ggplot2::ggplot(dsum, ggplot2::aes(x = t, y = median)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "#D55E00", alpha = 0.25
      ) +
      ggplot2::geom_line(colour = "#D55E00") +
      {if (!is.null(truth) && !is.null(truth$fprime))
        ggplot2::geom_line(
          data = data.frame(t = t_unique, value = truth$fprime),
          ggplot2::aes(x = t, y = value), inherit.aes = FALSE,
          linetype = 2
        )
      else NULL} +
      ggplot2::labs(title = titles[2L], x = "t", y = "f'(t)") +
      ggplot2::theme_minimal()
  }
  out
}
