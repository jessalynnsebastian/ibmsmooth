# Simulate integrated Brownian motion with interval-specific diffusion

Simulates the exact augmented-state transition used by both Stan models.

## Usage

``` r
simulate_ibm(t, lambda = 1, initial = c(0, 0), seed = NULL)
```

## Arguments

- t:

  Strictly increasing locations.

- lambda:

  Positive derivative diffusion variance rate, either scalar or one
  value per interval.

- initial:

  Initial state `(fprime, f)`.

- seed:

  Optional random seed.

## Value

A list containing `f`, `fprime`, `lambda_interval`, and the full
augmented state.
