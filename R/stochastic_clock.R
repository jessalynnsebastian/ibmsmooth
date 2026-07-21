#' Exact stochastic-clock IBM utilities
#'
#' These helpers expose the mathematics used by the Stan models. Operational
#' time over an observed interval is `q = lambda * dt` and the state is
#' `(B(T), I(T))`.
#'
#' @param q Positive operational-time increment.
#' @return A 2 by 2 matrix.
#' @export
ibm_clock_cholesky <- function(q) {
  if (length(q) != 1L || !is.finite(q) || q <= 0) stop("q must be positive")
  matrix(c(sqrt(q), q^(3 / 2) / 2, 0, q^(3 / 2) / sqrt(12)), 2L, 2L)
}

#' @rdname ibm_clock_cholesky
#' @export
ibm_clock_covariance <- function(q) {
  if (length(q) != 1L || !is.finite(q) || q <= 0) stop("q must be positive")
  matrix(c(q, q^2 / 2, q^2 / 2, q^3 / 3), 2L, 2L)
}

#' Simulate a literal time-changed integrated Brownian motion
#'
#' @param t Strictly increasing observed times.
#' @param lambda A positive clock rate, scalar or one value per interval.
#' @param initial Initial `(B, I)` state.
#' @param seed Optional random seed.
#' @return A list containing the operational state, clock quantities, function,
#' and left/right observed-time derivatives.
#' @export
simulate_ibm_clock <- function(t, lambda = 1, initial = c(0, 0), seed = NULL) {
  if (!is.numeric(t) || length(t) < 2L || any(!is.finite(t)) || any(diff(t) <= 0))
    stop("t must be a strictly increasing finite vector")
  n <- length(t); dt <- diff(t)
  if (length(lambda) == 1L) lambda <- rep(lambda, n - 1L)
  if (length(lambda) != n - 1L || any(!is.finite(lambda)) || any(lambda <= 0))
    stop("lambda must be positive and scalar or have length length(t) - 1")
  if (!is.null(seed)) set.seed(seed)
  q <- lambda * dt
  state <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c("B_operational", "f")))
  state[1L, ] <- initial
  for (i in 2:n) {
    mu <- c(state[i - 1L, 1L], state[i - 1L, 2L] + q[i - 1L] * state[i - 1L, 1L])
    state[i, ] <- mu + as.numeric(ibm_clock_cholesky(q[i - 1L]) %*% stats::rnorm(2L))
  }
  left <- c(NA_real_, lambda * state[-1L, 1L])
  right <- c(lambda * state[-n, 1L], NA_real_)
  list(t = t, state = state, B_operational = state[, 1L], f = state[, 2L],
       lambda_interval = lambda, q_interval = q,
       fprime_left = left, fprime_right = right)
}

#' Plot posterior stochastic-clock intervals
#'
#' @param ibmfit A stochastic-clock `ibmfit` object.
#' @param quantity One of interval clock rate, acceleration, or operational time.
#' @param probs Lower, center, and upper posterior probabilities.
#' @return A ggplot object using interval-step geometry.
#' @export
plot_clock <- function(ibmfit, quantity = c("lambda_interval", "r_interval", "q_interval"),
                       probs = c(0.025, 0.5, 0.975)) {
  quantity <- match.arg(quantity)
  if (is.null(ibmfit$stanfit)) stop("plot_clock currently requires a Stan fit")
  draws <- rstan::extract(ibmfit$stanfit, pars = quantity)[[quantity]]
  qs <- apply(draws, 2L, stats::quantile, probs = probs)
  grid <- ibmfit$data$time_grid_raw
  d <- data.frame(left = grid[-length(grid)], right = grid[-1L],
                  lower = qs[1L, ], center = qs[2L, ], upper = qs[3L, ])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = left, xend = right)) +
    ggplot2::geom_segment(ggplot2::aes(y = lower, yend = lower), alpha = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(y = upper, yend = upper), alpha = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(y = center, yend = center), linewidth = 0.8) +
    ggplot2::labs(x = "t", y = quantity) + ggplot2::theme_minimal()
  if (quantity == "lambda_interval" && ibmfit$adaptive == "baseline_horseshoe") {
    base <- rstan::extract(ibmfit$stanfit, pars = "lambda0")$lambda0
    p <- p + ggplot2::geom_hline(yintercept = stats::median(base), linetype = 2)
  }
  p
}
