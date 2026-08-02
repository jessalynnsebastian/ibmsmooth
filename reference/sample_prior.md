# Draw from the prior associated with a fitted model

The prior fit uses the same locations, scaling, and prior settings as
`ibmfit`, but omits the observation likelihood.

## Usage

``` r
sample_prior(
  ibmfit,
  iter = 1000,
  chains = 2,
  cores = getOption("mc.cores", chains),
  ...
)
```

## Arguments

- ibmfit:

  An `ibmfit` object returned by
  [`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md).

- iter, chains, cores:

  Stan sampling controls for the prior fit.

- ...:

  Additional arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html),
  such as `seed` or `refresh`.

## Value

An `ibmfit` containing prior draws.
