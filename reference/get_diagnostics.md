# Extract computation and MCMC diagnostics

Reports elapsed sampling time, effective sample sizes, split R-hat,
divergences, leapfrog counts, treedepth, acceptance statistics, and
energy Bayesian fraction of missing information (E-BFMI).

## Usage

``` r
get_diagnostics(ibmfit, pars = NULL)
```

## Arguments

- ibmfit:

  An `ibmfit`.

- pars:

  Optional character vector of Stan parameters. By default all monitored
  parameters are included.

## Value

A list with `timing`, `sampler`, `parameters`, and `overview`
components.
