data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  // OBSERVATION_DATA
  real<lower=0> global_scale;
  int<lower=1, upper=2> global_prior;
  real<lower=0> initial_sd;
  int<lower=0, upper=1> prior_only;
}
parameters {
  // OBSERVATION_PARAMETERS
  real<lower=0, upper=1> gamma_unif;
  vector<lower=0, upper=1>[T - 1] xi_unif;
  vector[T] fprime;
  vector[T] f;
}
transformed parameters {
  // OBSERVATION_TRANSFORMED_PARAMETERS
  real<lower=0> gamma = global_prior == 1
    ? global_scale * tan(0.5 * pi() * gamma_unif)
    : global_scale * inv_Phi(0.5 + 0.5 * gamma_unif);
  vector<lower=0>[T - 1] xi = tan(0.5 * pi() * xi_unif);
}
model {
  // OBSERVATION_PRIORS
  gamma_unif ~ uniform(0, 1);
  xi_unif ~ uniform(0, 1);
  fprime[1] ~ normal(0, initial_sd);
  f[1] ~ normal(0, initial_sd);

  for (i in 2:T) {
    real h = deltat[i - 1];
    real tau_i = gamma * xi[i - 1];
    vector[2] state = [fprime[i], f[i]]';
    vector[2] transition_mean =
      [fprime[i - 1], f[i - 1] + h * fprime[i - 1]]';
    matrix[2, 2] transition_cholesky;

    transition_cholesky[1, 1] = tau_i * sqrt(h);
    transition_cholesky[1, 2] = 0;
    transition_cholesky[2, 1] = tau_i * h * sqrt(h) / 2;
    transition_cholesky[2, 2] =
      tau_i * h * sqrt(h) * inv_sqrt(12.0);
    state ~ multi_normal_cholesky(transition_mean, transition_cholesky);
  }
  if (!prior_only) {
    // OBSERVATION_LIKELIHOOD
  }
}
generated quantities {
  vector<lower=0>[T - 1] lambda_interval = square(gamma * xi);
}
