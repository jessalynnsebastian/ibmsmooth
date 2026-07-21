# Package index

## Fit models

- [`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md)
  : Unified interface for integrated Brownian motion smoothing
- [`ibm_smooth()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_smooth.md)
  : Fit an integrated Brownian motion (IBM) smoother using Stan
- [`ibm_inla_fit()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_inla_fit.md)
  : Fit an integrated Brownian motion smoother using INLA
- [`ibm_horseshoe_marginalized()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_horseshoe_marginalized.md)
  : Fit a marginalized Gaussian horseshoe IBM
- [`validate_marginalized_horseshoe()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/validate_marginalized_horseshoe.md)
  : Validate the marginalized Kalman likelihood against a dense Gaussian
  result

## Summarize posterior draws

- [`get_curve_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_curve_samples.md)
  : Extract posterior curve samples
- [`get_curve_summary()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_curve_summary.md)
  : Summarize posterior curves
- [`get_hyperparameter_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_hyperparameter_samples.md)
  : Extract posterior hyperparameter samples
- [`get_hyperparameter_summary()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_hyperparameter_summary.md)
  : Summarize posterior hyperparameters
- [`get_samples()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_samples.md)
  : Extract posterior samples from an IBM fit
- [`get_samples_inla()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_samples_inla.md)
  : Extract posterior samples from an INLA fit

## Plot and print results

- [`plot(`*`<ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot.ibmfit.md)
  : Plot method for ibmfit objects
- [`plot_curve()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_curve.md)
  : Plot posterior function and derivative summaries
- [`plot_precision_curve()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_precision_curve.md)
  : Plot the locally adaptive process precision curve
- [`summary(`*`<ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/summary.ibmfit.md)
  : Summarize an ibmfit object
- [`print(`*`<summary_ibmfit>`*`)`](https://jessalynnsebastian.github.io/ibmsmooth/reference/print.summary_ibmfit.md)
  : Print an ibmfit summary
