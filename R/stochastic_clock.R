#' Integrated Brownian motion transition matrices
#'
#' @param delta Positive interval length.
#' @param lambda Positive derivative diffusion variance rate.
#' @return A 2 by 2 covariance matrix or its lower Cholesky factor, ordered as
#'   `(fprime, f)`.
#' @export
ibm_transition_covariance <- function(delta, lambda = 1) {
  .check_transition_args(delta, lambda)
  lambda * matrix(
    c(delta, delta^2 / 2, delta^2 / 2, delta^3 / 3),
    2L, 2L
  )
}

#' @rdname ibm_transition_covariance
#' @export
ibm_transition_cholesky <- function(delta, lambda = 1) {
  .check_transition_args(delta, lambda)
  sqrt(lambda) * matrix(
    c(sqrt(delta), delta^(3 / 2) / 2,
      0, delta^(3 / 2) / sqrt(12)),
    2L, 2L
  )
}

.check_transition_args <- function(delta, lambda) {
  if (!is.numeric(delta) || length(delta) != 1L ||
      !is.finite(delta) || delta <= 0 ||
      !is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda <= 0) {
    stop("delta and lambda must be positive finite scalars.", call. = FALSE)
  }
}

#' Simulate integrated Brownian motion with interval-specific diffusion
#'
#' Simulates the exact augmented-state transition used by both Stan models.
#'
#' @param t Strictly increasing locations.
#' @param lambda Positive derivative diffusion variance rate, either scalar or
#'   one value per interval.
#' @param initial Initial state `(fprime, f)`.
#' @param seed Optional random seed.
#' @return A list containing `f`, `fprime`, `lambda_interval`, and the full
#'   augmented state.
#' @export
simulate_ibm <- function(t, lambda = 1, initial = c(0, 0), seed = NULL) {
  if (!is.numeric(t) || length(t) < 2L || any(!is.finite(t)) ||
      any(diff(t) <= 0)) {
    stop("t must be a strictly increasing finite numeric vector.", call. = FALSE)
  }
  n <- length(t)
  if (length(lambda) == 1L) lambda <- rep(lambda, n - 1L)
  if (!is.numeric(lambda) || length(lambda) != n - 1L ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must be positive and scalar or have length length(t) - 1.",
         call. = FALSE)
  }
  if (!is.numeric(initial) || length(initial) != 2L ||
      any(!is.finite(initial))) {
    stop("initial must be a finite numeric vector of length two.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  state <- matrix(NA_real_, n, 2L,
                  dimnames = list(NULL, c("fprime", "f")))
  state[1L, ] <- initial
  for (i in 2:n) {
    h <- t[i] - t[i - 1L]
    mean <- c(state[i - 1L, 1L],
              state[i - 1L, 2L] + h * state[i - 1L, 1L])
    state[i, ] <- mean +
      as.numeric(ibm_transition_cholesky(h, lambda[i - 1L]) %*%
                   stats::rnorm(2L))
  }
  list(
    t = t, state = state, fprime = state[, 1L], f = state[, 2L],
    lambda_interval = lambda
  )
}

#' Plot interval-specific diffusion rates
#'
#' @param ibmfit A locally adaptive `ibmfit`.
#' @param probs Lower, center, and upper posterior probabilities.
#' @return A ggplot object.
#' @export
plot_diffusion <- function(ibmfit, probs = c(0.025, 0.5, 0.975)) {
  if (!inherits(ibmfit, "ibmfit") || !isTRUE(ibmfit$adaptive)) {
    stop("ibmfit must be a locally adaptive fit.", call. = FALSE)
  }
  if (length(probs) != 3L || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || any(diff(probs) <= 0)) {
    stop("probs must contain increasing lower, center, and upper probabilities.",
         call. = FALSE)
  }
  draws <- get_samples(ibmfit, "lambda_interval")
  qs <- apply(draws, 2L, stats::quantile, probs = probs)
  grid <- ibmfit$data$time_grid_raw
  d <- data.frame(
    left = grid[-length(grid)], right = grid[-1L],
    lower = qs[1L, ], center = qs[2L, ], upper = qs[3L, ]
  )
  ggplot2::ggplot(d) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = left, xmax = right, ymin = lower, ymax = upper),
      fill = "#56B4E9", alpha = 0.25
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = left, xend = right, y = center, yend = center),
      colour = "#0072B2", linewidth = 0.9
    ) +
    ggplot2::labs(
      title = "Posterior interval diffusion variance rate",
      x = "t", y = expression(lambda[j])
    ) +
    ggplot2::theme_minimal()
}
