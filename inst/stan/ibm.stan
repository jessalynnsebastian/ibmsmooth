data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  vector[N_obs] y_obs;
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;
  real<lower=0> initial_sd;
}
parameters {
  real log_sigma_raw;
  real<lower=0> tau;
  vector[2] z_initial;
  array[T - 1] vector[2] z_transition;
}
transformed parameters {
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd * log_sigma_raw);
  vector[T] fprime;
  vector[T] f;

  fprime[1] = initial_sd * z_initial[1];
  f[1] = initial_sd * z_initial[2];
  for (i in 2:T) {
    real h = deltat[i - 1];
    real h32 = h * sqrt(h);
    fprime[i] = fprime[i - 1]
      + tau * sqrt(h) * z_transition[i - 1][1];
    f[i] = f[i - 1] + h * fprime[i - 1]
      + tau * h32
        * (0.5 * z_transition[i - 1][1]
           + inv_sqrt(12.0) * z_transition[i - 1][2]);
  }
}
model {
  log_sigma_raw ~ std_normal();
  z_initial ~ std_normal();
  for (i in 1:(T - 1)) z_transition[i] ~ std_normal();
  // SMOOTHING_PRIOR
  for (n in 1:N_obs)
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
}
