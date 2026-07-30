#' Check whether a regularized horseshoe retry may help
#'
#' Looks for the specific combination of poor HMC diagnostics and one or more
#' local half-Cauchy coordinates concentrated near their upper boundary. This
#' is the failure mode in which a finite slab is most likely to help.
#'
#' @param ibmfit An adaptive ordinary-horseshoe `ibmfit`.
#' @param boundary Median `xi_unif` threshold used to identify an extreme
#'   local scale.
#' @param rhat_threshold R-hat threshold.
#' @param ess_threshold Effective sample size threshold.
#' @return A list containing `recommended`, `reasons`, sampler diagnostics, and
#'   the indices and posterior medians of extreme local scales.
#' @export
check_regularized_horseshoe <- function(
    ibmfit, boundary = 0.99, rhat_threshold = 1.01,
    ess_threshold = 100) {
  .check_ibmfit(ibmfit)
  if (!isTRUE(ibmfit$adaptive) ||
      !identical(ibmfit$adaptive_prior, "horseshoe")) {
    return(list(
      recommended = FALSE,
      reasons = "The fit is not an adaptive ordinary-horseshoe fit.",
      extreme_intervals = integer(),
      extreme_medians = numeric()
    ))
  }
  for (x in list(boundary, rhat_threshold, ess_threshold)) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop("Diagnostic thresholds must be finite numeric scalars.",
           call. = FALSE)
    }
  }
  if (boundary <= 0 || boundary >= 1 ||
      rhat_threshold <= 1 || ess_threshold <= 0) {
    stop("Diagnostic thresholds are outside their valid ranges.",
         call. = FALSE)
  }

  fit <- ibmfit$stanfit
  xi_unif <- rstan::extract(
    fit, pars = "xi_unif", permuted = TRUE
  )$xi_unif
  if (is.null(dim(xi_unif))) xi_unif <- matrix(xi_unif, ncol = 1L)
  local_medians <- apply(xi_unif, 2L, stats::median)
  extreme <- which(local_medians >= boundary)

  parameter_summary <- rstan::summary(
    fit, pars = c("gamma_unif", "xi_unif", "z_transition", "f", "fprime")
  )$summary
  max_rhat <- suppressWarnings(max(
    parameter_summary[, "Rhat"], na.rm = TRUE
  ))
  min_ess <- suppressWarnings(min(
    parameter_summary[, "n_eff"], na.rm = TRUE
  ))
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler, function(x) sum(x[, "divergent__"]), numeric(1L)
  ))
  max_depth <- ibmfit$fit_spec$max_treedepth
  treedepth_hits <- if (length(max_depth)) {
    sum(vapply(
      sampler,
      function(x) sum(x[, "treedepth__"] >= max_depth),
      numeric(1L)
    ))
  } else {
    NA_integer_
  }
  poor_sampling <- divergences > 0L ||
    (!is.na(treedepth_hits) && treedepth_hits > 0L) ||
    is.finite(max_rhat) && max_rhat > rhat_threshold ||
    is.finite(min_ess) && min_ess < ess_threshold

  reasons <- character()
  if (length(extreme)) {
    reasons <- c(
      reasons,
      paste0(
        "median xi_unif exceeded ", boundary, " at interval",
        if (length(extreme) > 1L) "s " else " ",
        paste(extreme, collapse = ", ")
      )
    )
  }
  if (divergences > 0L) {
    reasons <- c(reasons, paste(divergences, "divergent transitions"))
  }
  if (!is.na(treedepth_hits) && treedepth_hits > 0L) {
    reasons <- c(reasons, paste(treedepth_hits, "maximum-treedepth hits"))
  }
  if (is.finite(max_rhat) && max_rhat > rhat_threshold) {
    reasons <- c(reasons, paste0("maximum R-hat was ", signif(max_rhat, 4)))
  }
  if (is.finite(min_ess) && min_ess < ess_threshold) {
    reasons <- c(reasons, paste0("minimum ESS was ", signif(min_ess, 4)))
  }

  list(
    recommended = length(extreme) > 0L && isTRUE(poor_sampling),
    reasons = reasons,
    extreme_intervals = extreme,
    extreme_medians = local_medians[extreme],
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    max_rhat = max_rhat,
    min_ess = min_ess,
    thresholds = list(
      boundary = boundary, rhat = rhat_threshold, ess = ess_threshold
    )
  )
}
