data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  vector[N_obs] y_obs;
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;
  real<lower=0> initial_sd;
  int<lower=0, upper=1> prior_only;
}
parameters {
  real log_sigma_raw;
  real<lower=0> tau;
  vector[T] fprime;
  vector[T] f;
}
transformed parameters {
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd * log_sigma_raw);
}
model {
  log_sigma_raw ~ std_normal();
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
  if (!prior_only)
    for (n in 1:N_obs)
      y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
}
