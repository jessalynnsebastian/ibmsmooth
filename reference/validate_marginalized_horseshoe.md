# Validate the marginalized Kalman likelihood against a dense Gaussian result

Intended as a small-data numerical check for the marginalized model.

## Usage

``` r
validate_marginalized_horseshoe(t, y, sigma, tau, initial_sd = 5)
```

## Arguments

- t, y:

  Observation times and scaled observations. Times must be unique.

- sigma:

  Observation standard deviation.

- tau:

  Vector of transition scales of length `length(t) - 1`.

- initial_sd:

  Scale of the correlated unit-time IBM initial-state prior.

## Value

The Kalman and dense log likelihoods and their difference.
