#' Extract computation and MCMC diagnostics
#'
#' Reports elapsed sampling time, effective sample sizes, split R-hat,
#' divergences, maximum-treedepth hits, acceptance statistics, and
#' energy Bayesian fraction of missing information (E-BFMI).
#'
#' @param ibmfit An `ibmfit`.
#' @param pars Optional character vector of Stan parameters. By default all
#'   monitored parameters are included.
#' @return A list with `timing`, `sampler`, `parameters`, and `overview`
#'   components.
#' @export
get_diagnostics <- function(ibmfit, pars = NULL) {
  .check_ibmfit(ibmfit)
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required to extract diagnostics.",
         call. = FALSE)
  }
  fit <- ibmfit$stanfit
  elapsed <- rstan::get_elapsed_time(fit)
  timing <- data.frame(
    chain = seq_len(nrow(elapsed)),
    warmup_seconds = elapsed[, "warmup"],
    sampling_seconds = elapsed[, "sample"],
    total_seconds = rowSums(elapsed),
    row.names = NULL
  )
  timing$cpu_seconds_total <- sum(timing$total_seconds)
  timing$wall_seconds_approx <- max(timing$total_seconds)

  sampler_raw <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  max_depth <- ibmfit$fit_spec$max_treedepth
  chain_diagnostics <- lapply(seq_along(sampler_raw), function(i) {
    x <- sampler_raw[[i]]
    energy <- x[, "energy__"]
    ebfmi <- if (length(energy) > 1L && stats::var(energy) > 0) {
      mean(diff(energy)^2) / stats::var(energy)
    } else {
      NA_real_
    }
    data.frame(
      chain = i,
      iterations = nrow(x),
      divergences = sum(x[, "divergent__"]),
      treedepth_hits = if (length(max_depth)) {
        sum(x[, "treedepth__"] >= max_depth)
      } else NA_integer_,
      mean_accept_stat = mean(x[, "accept_stat__"]),
      mean_stepsize = mean(x[, "stepsize__"]),
      ebfmi = ebfmi
    )
  })
  sampler <- do.call(rbind, chain_diagnostics)

  sm <- rstan::summary(fit, pars = pars)$summary
  parameters <- data.frame(
    parameter = rownames(sm), sm, row.names = NULL, check.names = FALSE
  )
  overview <- data.frame(
    cpu_seconds = sum(timing$total_seconds),
    wall_seconds_approx = max(timing$total_seconds),
    min_ess = suppressWarnings(min(parameters$n_eff, na.rm = TRUE)),
    max_rhat = suppressWarnings(max(parameters$Rhat, na.rm = TRUE)),
    divergences = sum(sampler$divergences),
    treedepth_hits = sum(sampler$treedepth_hits, na.rm = TRUE),
    min_ebfmi = suppressWarnings(min(sampler$ebfmi, na.rm = TRUE))
  )
  list(timing = timing, sampler = sampler, parameters = parameters,
       overview = overview)
}
