#' Unified interface for integrated Brownian motion smoothing
#'
#' Fit an integrated Brownian motion smoother using either Stan or INLA.
#'
#' @param t Numeric vector of time points.
#' @param y Numeric vector of observations at times `t`.
#' @param infer_at Optional numeric vector of additional time points at which
#'   inference should be returned.
#' @param method Computational backend, either `"stan"` or `"inla"`.
#' @param adaptive Logical. If `FALSE`, use a single smoothness parameter. If
#'   `TRUE`, fit a locally adaptive model.
#' @param stan_adaptive_method Stan-only adaptive prior. The default is
#'   `"horseshoe"`.
#' @param stan_horseshoe_engine Stan engine used for the adaptive horseshoe
#'   model. `"joint"` retains the existing sampler; `"marginalized"` uses the
#'   Gaussian-data Kalman marginalization and a simulation smoother.
#' @param ... Additional arguments passed to `ibm_smooth()` or `ibm_inla_fit()`.
#'
#' @return An object of class `ibmfit`.
#' @export
ibm <- function(t, y, infer_at = NULL, method = c("stan", "inla"),
                adaptive = FALSE,
                stan_adaptive_method = c("horseshoe", "rw", "rhs", "bridge"),
                stan_horseshoe_engine = c("joint", "marginalized"),
                ...) {
  method <- match.arg(method)

  if (method == "stan") {
    stan_adaptive_method <- match.arg(stan_adaptive_method)
    stan_horseshoe_engine <- match.arg(stan_horseshoe_engine)
    if (isTRUE(adaptive) && stan_adaptive_method == "horseshoe" &&
        stan_horseshoe_engine == "marginalized") {
      return(ibm_horseshoe_marginalized(
        t = t, y = y, infer_at = infer_at, ...
      ))
    }
    stan_adapt <- if (isTRUE(adaptive)) stan_adaptive_method else "nonadaptive"
    return(ibm_smooth(
      t = t,
      y = y,
      infer_at = infer_at,
      adaptive = stan_adapt,
      ...
    ))
  }

  ibm_inla_fit(
    t = t,
    y = y,
    infer_at = infer_at,
    adaptive = adaptive,
    ...
  )
}
