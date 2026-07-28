.ibm_A <- function(h) matrix(c(1, h, 0, 1), 2L, 2L)

.ibm_Q <- function(h) {
  matrix(c(h, h^2 / 2, h^2 / 2, h^3 / 3), 2L, 2L)
}

.draw_bivariate <- function(mean, covariance) {
  covariance <- (covariance + t(covariance)) / 2
  eigen_decomp <- eigen(covariance, symmetric = TRUE)
  values <- pmax(eigen_decomp$values, 0)
  mean + as.numeric(
    eigen_decomp$vectors %*% (sqrt(values) * stats::rnorm(2L))
  )
}

#' Predict an IBM curve after fitting
#'
#' Draws function values and derivatives at arbitrary locations without adding
#' those locations to the Stan latent state. Within the observed range, draws
#' use the exact IBM bridge conditional on adjacent sampled states. Outside the
#' range they use forward or reverse IBM transitions.
#'
#' @param ibmfit An `ibmfit`.
#' @param new_t Finite numeric prediction locations. The default uses the
#'   observed locations together with the fit's original `infer_at`.
#' @param n_samples Maximum number of posterior draws.
#' @param seed Optional random seed for bridge and extrapolation draws.
#' @return A list with `t`, `f`, and `fprime`. Draws are rows and locations are
#'   columns.
#' @export
predict_curve <- function(ibmfit, new_t = NULL, n_samples = 1000,
                          seed = NULL) {
  .check_ibmfit(ibmfit)
  latent_t <- ibmfit$data$time_grid_raw
  if (is.null(new_t)) {
    new_t <- ibmfit$data$prediction_grid_raw
    if (is.null(new_t)) new_t <- latent_t
  }
  if (!is.numeric(new_t) || !length(new_t) || any(!is.finite(new_t))) {
    stop("new_t must be a non-empty finite numeric vector.", call. = FALSE)
  }
  new_t <- sort(unique(as.numeric(new_t)))
  f <- get_samples(ibmfit, "f", n_samples)
  fprime <- get_samples(ibmfit, "fprime", n_samples)
  n_draws <- nrow(f)
  if (!is.null(seed)) set.seed(seed)

  result_f <- result_fprime <- matrix(
    NA_real_, nrow = n_draws, ncol = length(new_t)
  )
  exact <- match(new_t, latent_t)
  exact_columns <- which(!is.na(exact))
  if (length(exact_columns)) {
    result_f[, exact_columns] <- f[, exact[exact_columns], drop = FALSE]
    result_fprime[, exact_columns] <-
      fprime[, exact[exact_columns], drop = FALSE]
  }
  if (length(exact_columns) == length(new_t)) {
    colnames(result_f) <- colnames(result_fprime) <- as.character(new_t)
    return(list(t = new_t, f = result_f, fprime = result_fprime))
  }

  if (isTRUE(ibmfit$adaptive)) {
    lambda <- get_samples(ibmfit, "lambda_interval", n_draws)
  } else {
    tau <- get_samples(ibmfit, "tau", n_draws)[, 1L]
    lambda <- matrix(tau^2, nrow = n_draws,
                     ncol = length(latent_t) - 1L)
  }

  for (draw in seq_len(n_draws)) {
    states <- rbind(fprime[draw, ], f[draw, ])

    # Reverse-time extrapolation to the left of the first observed state.
    left_queries <- sort(new_t[new_t < latent_t[1L]], decreasing = TRUE)
    current <- states[, 1L]
    current_t <- latent_t[1L]
    for (query in left_queries) {
      h <- current_t - query
      reverse_A <- .ibm_A(-h)
      covariance <- lambda[draw, 1L] *
        reverse_A %*% .ibm_Q(h) %*% t(reverse_A)
      current <- .draw_bivariate(reverse_A %*% current, covariance)
      column <- match(query, new_t)
      result_fprime[draw, column] <- current[1L]
      result_f[draw, column] <- current[2L]
      current_t <- query
    }

    # Conditional bridges within each observed interval.
    for (interval in seq_len(length(latent_t) - 1L)) {
      left_t <- latent_t[interval]
      right_t <- latent_t[interval + 1L]
      queries <- new_t[new_t > left_t & new_t < right_t]
      if (!length(queries)) next
      current <- states[, interval]
      current_t <- left_t
      right_state <- states[, interval + 1L]
      for (query in queries) {
        a <- query - current_t
        b <- right_t - query
        total <- a + b
        A_a <- .ibm_A(a)
        A_b <- .ibm_A(b)
        A_total <- .ibm_A(total)
        cross_covariance <- .ibm_Q(a) %*% t(A_b)
        gain <- cross_covariance %*% solve(.ibm_Q(total))
        conditional_mean <- A_a %*% current +
          gain %*% (right_state - A_total %*% current)
        conditional_covariance <- lambda[draw, interval] *
          (.ibm_Q(a) - gain %*% t(cross_covariance))
        current <- .draw_bivariate(
          conditional_mean, conditional_covariance
        )
        column <- match(query, new_t)
        result_fprime[draw, column] <- current[1L]
        result_f[draw, column] <- current[2L]
        current_t <- query
      }
    }

    # Forward extrapolation to the right of the last observed state.
    right_queries <- new_t[new_t > latent_t[length(latent_t)]]
    current <- states[, length(latent_t)]
    current_t <- latent_t[length(latent_t)]
    for (query in right_queries) {
      h <- query - current_t
      current <- .draw_bivariate(
        .ibm_A(h) %*% current,
        lambda[draw, ncol(lambda)] * .ibm_Q(h)
      )
      column <- match(query, new_t)
      result_fprime[draw, column] <- current[1L]
      result_f[draw, column] <- current[2L]
      current_t <- query
    }
  }
  colnames(result_f) <- colnames(result_fprime) <- as.character(new_t)
  list(t = new_t, f = result_f, fprime = result_fprime)
}
