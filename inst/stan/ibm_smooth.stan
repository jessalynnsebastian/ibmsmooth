data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T-1] deltat;
  vector[N_obs] y_obs;
  int<lower=0, upper=1> regular; // retained for API compatibility
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;
  real log_tau_mu;               // prior for log clock rate (legacy data name)
  real<lower=0> log_tau_sd;
  real<lower=0> initial_sd;
}
parameters {
  real log_sigma_raw;
  real log_lambda0_raw;
  matrix[2, T] z_state;
}
transformed parameters {
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd * log_sigma_raw);
  real<lower=0> lambda0 = exp(log_tau_mu + log_tau_sd * log_lambda0_raw);
  vector[T-1] lambda_interval = rep_vector(lambda0, T-1);
  vector[T-1] q_interval = lambda_interval .* deltat;
  vector[T] B_operational;
  vector[T] f;
  vector[T] fprime_left;
  vector[T] fprime_right;

  B_operational[1] = initial_sd * z_state[1, 1];
  f[1] = initial_sd * z_state[2, 1];
  for (i in 2:T) {
    real q = q_interval[i-1];
    real sqrt_q = sqrt(q);
    real q32 = q * sqrt_q;
    B_operational[i] = B_operational[i-1] + sqrt_q * z_state[1, i];
    f[i] = f[i-1] + q * B_operational[i-1]
      + 0.5 * q32 * z_state[1, i]
      + q32 * inv_sqrt(12.0) * z_state[2, i];
  }
  // Boundary placeholders are masked to NA by the R extractor.
  fprime_left[1] = 0;
  fprime_right[T] = 0;
  for (i in 2:T) fprime_left[i] = lambda_interval[i-1] * B_operational[i];
  for (i in 1:(T-1)) fprime_right[i] = lambda_interval[i] * B_operational[i];
}
model {
  log_sigma_raw ~ std_normal();
  log_lambda0_raw ~ std_normal();
  to_vector(z_state) ~ std_normal();
  for (n in 1:N_obs) y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
}
