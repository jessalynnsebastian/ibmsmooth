# Generate posterior predictive observations

Generate posterior predictive observations

## Usage

``` r
posterior_predict(ibmfit, n_samples = 500, seed = NULL)
```

## Arguments

- ibmfit:

  A posterior `ibmfit`.

- n_samples:

  Number of posterior predictive data sets.

- seed:

  Optional random seed.

## Value

A matrix with predictive data sets in rows and observations in columns,
in the original response units.
