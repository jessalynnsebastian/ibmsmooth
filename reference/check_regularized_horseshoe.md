# Check whether a regularized horseshoe retry may help

Looks for the specific combination of poor HMC diagnostics and one or
more local half-Cauchy coordinates concentrated near their upper
boundary. This is the failure mode in which a finite slab is most likely
to help.

## Usage

``` r
check_regularized_horseshoe(
  ibmfit,
  boundary = 0.99,
  rhat_threshold = 1.01,
  ess_threshold = 100
)
```

## Arguments

- ibmfit:

  An adaptive ordinary-horseshoe `ibmfit`.

- boundary:

  Median `xi_unif` threshold used to identify an extreme local scale.

- rhat_threshold:

  R-hat threshold.

- ess_threshold:

  Effective sample size threshold.

## Value

A list containing `recommended`, `reasons`, sampler diagnostics, and the
indices and posterior medians of extreme local scales.
