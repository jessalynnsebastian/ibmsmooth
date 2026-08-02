# Package index

## Fit models

- [`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md)
  : Fit an integrated Brownian motion smoother
- [`ibm_smooth()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_smooth.md)
  : Fit an IBM smoother using Stan
- [`sample_prior()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/sample_prior.md)
  : Draw from the prior associated with a fitted model

## Prediction

- [`predict_curve()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/predict_curve.md)
  : Predict an IBM curve after fitting
- [`posterior_predict()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/posterior_predict.md)
  : Generate posterior predictive observations

## Exact IBM transitions and simulation

- [`ibm_transition_covariance()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_transition_covariance.md)
  [`ibm_transition_cholesky()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_transition_covariance.md)
  : Integrated Brownian motion transition matrices
- [`simulate_ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/simulate_ibm.md)
  : Simulate integrated Brownian motion with interval-specific diffusion

## Model diagnostics

- [`check_regularized_horseshoe()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/check_regularized_horseshoe.md)
  : Check whether a regularized horseshoe retry may help

## Summarize and plot

- [`summarize_performance()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/summarize_performance.md)
  : Summarize function and derivative estimation performance
- [`get_curve_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_curve_samples.md)
  : Extract posterior curve samples
- [`get_curve_summary()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_curve_summary.md)
  : Summarize posterior curves
- [`get_diagnostics()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_diagnostics.md)
  : Extract computation and MCMC diagnostics
- [`get_hyperparameter_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_hyperparameter_samples.md)
  : Extract posterior hyperparameter samples
- [`get_hyperparameter_summary()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_hyperparameter_summary.md)
  : Summarize posterior hyperparameters
- [`get_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_samples.md)
  : Extract posterior samples
- [`plot(`*`<ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot.ibmfit.md)
  : Plot an integrated Brownian motion fit
- [`plot_curve()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_curve.md)
  : Plot posterior function and derivative summaries
- [`plot_diffusion()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_diffusion.md)
  : Plot interval-specific diffusion rates
- [`plot_pp_check()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_pp_check.md)
  : Visual posterior predictive check
- [`plot_prior_posterior()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_prior_posterior.md)
  : Plot prior and posterior distributions
- [`summary(`*`<ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/summary.ibmfit.md)
  : Summarize an IBM fit
- [`print(`*`<summary_ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/print.summary_ibmfit.md)
  : Print an IBM fit summary
