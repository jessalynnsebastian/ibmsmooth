#' Verify the standard and fast IBM prior representations
#'
#' Compares the joint covariance of function values obtained from the
#' augmented two-dimensional state transition with that obtained after
#' analytically filtering out the derivative. Because both representations are
#' zero-mean Gaussian conditional on `lambda`, matching this covariance verifies
#' equality of their complete conditional distributions.
#'
#' @param t Strictly increasing finite locations.
#' @param lambda Positive derivative diffusion variance rate, either scalar or
#'   one value per interval.
#' @param initial_sd Positive standard deviation for both initial state
#'   components.
#' @param tolerance Numerical tolerance used for the `equivalent` result.
#' @return A list containing `equivalent`, `max_abs_covariance_difference`, and
#'   the covariance matrices from both representations.
#' @export
verify_ibm_equivalence <- function(t, lambda = 1, initial_sd = 5,
                                   tolerance = 1e-10) {
  if (!is.numeric(t) || length(t) < 2L || any(!is.finite(t)) ||
      any(diff(t) <= 0)) {
    stop("t must be a strictly increasing finite numeric vector.",
         call. = FALSE)
  }
  intervals <- length(t) - 1L
  if (length(lambda) == 1L) lambda <- rep(lambda, intervals)
  if (!is.numeric(lambda) || length(lambda) != intervals ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must be positive and scalar or have length length(t) - 1.",
         call. = FALSE)
  }
  if (!is.numeric(initial_sd) || length(initial_sd) != 1L ||
      !is.finite(initial_sd) || initial_sd <= 0) {
    stop("initial_sd must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive finite scalar.", call. = FALSE)
  }

  n <- length(t)
  augmented <- matrix(0, 2L * n, 2L * n)
  augmented[1:2, 1:2] <- diag(initial_sd^2, 2L)
  for (j in 2:n) {
    h <- t[j] - t[j - 1L]
    transition <- matrix(c(1, h, 0, 1), 2L, 2L)
    previous <- (2L * j - 3L):(2L * j - 2L)
    current <- (2L * j - 1L):(2L * j)
    augmented[current, current] <-
      transition %*% augmented[previous, previous] %*% t(transition) +
      ibm_transition_covariance(h, lambda[j - 1L])
    earlier <- seq_len(2L * j - 2L)
    augmented[current, earlier] <-
      transition %*% augmented[previous, earlier, drop = FALSE]
    augmented[earlier, current] <- t(augmented[current, earlier])
  }
  function_indices <- seq.int(2L, 2L * n, by = 2L)
  augmented_covariance <-
    augmented[function_indices, function_indices, drop = FALSE]

  filtered_draw <- function(z) {
    f <- numeric(n)
    f[1L] <- initial_sd * z[1L]
    derivative_mean <- 0
    derivative_variance <- initial_sd^2
    for (j in 2:n) {
      h <- t[j] - t[j - 1L]
      prediction_variance <- h^2 * derivative_variance +
        lambda[j - 1L] * h^3 / 3
      cross_covariance <- h * derivative_variance +
        lambda[j - 1L] * h^2 / 2
      innovation <- sqrt(prediction_variance) * z[j]
      f[j] <- f[j - 1L] + h * derivative_mean + innovation
      derivative_mean <- derivative_mean +
        cross_covariance / prediction_variance * innovation
      derivative_variance <- derivative_variance +
        lambda[j - 1L] * h -
        cross_covariance^2 / prediction_variance
    }
    f
  }
  innovation_map <- vapply(seq_len(n), function(j) {
    filtered_draw(replace(numeric(n), j, 1))
  }, numeric(n))
  filtered_covariance <- innovation_map %*% t(innovation_map)
  difference <- max(abs(augmented_covariance - filtered_covariance))

  list(
    equivalent = isTRUE(difference <= tolerance),
    max_abs_covariance_difference = difference,
    tolerance = tolerance,
    augmented_covariance = augmented_covariance,
    filtered_covariance = filtered_covariance
  )
}
