data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  // OBSERVATION_DATA
  real<lower=0> initial_sd;
  int<lower=0, upper=1> prior_only;
}
parameters {
  // OBSERVATION_PARAMETERS
  real<lower=0> tau;
  vector[T] fprime;
  vector[T] f;
}
transformed parameters {
  // OBSERVATION_TRANSFORMED_PARAMETERS
}
model {
  // OBSERVATION_PRIORS
  fprime[1] ~ normal(0, initial_sd);
  f[1] ~ normal(0, initial_sd);
  // SMOOTHING_PRIOR

  for (i in 2:T) {
    real h = deltat[i - 1];
    vector[2] state = [fprime[i], f[i]]';
    vector[2] transition_mean =
      [fprime[i - 1], f[i - 1] + h * fprime[i - 1]]';
    matrix[2, 2] transition_cholesky;

    transition_cholesky[1, 1] = tau * sqrt(h);
    transition_cholesky[1, 2] = 0;
    transition_cholesky[2, 1] = tau * h * sqrt(h) / 2;
    transition_cholesky[2, 2] = tau * h * sqrt(h) * inv_sqrt(12.0);
    state ~ multi_normal_cholesky(transition_mean, transition_cholesky);
  }
  if (!prior_only) {
    // OBSERVATION_LIKELIHOOD
  }
}
