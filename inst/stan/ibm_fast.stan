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
  real z_f_initial;
  vector[T - 1] z_f_transition;
}
transformed parameters {
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd * log_sigma_raw);
  vector[T] f;
  vector[T] derivative_filter_mean;
  vector<lower=0>[T] derivative_filter_variance;

  f[1] = initial_sd * z_f_initial;
  derivative_filter_mean[1] = 0;
  derivative_filter_variance[1] = square(initial_sd);
  for (i in 2:T) {
    real h = deltat[i - 1];
    real lambda_i = square(tau);
    real derivative_variance = derivative_filter_variance[i - 1];
    real prediction_variance =
      square(h) * derivative_variance + lambda_i * h * square(h) / 3;
    real cross_covariance =
      h * derivative_variance + lambda_i * square(h) / 2;
    real innovation = sqrt(prediction_variance) * z_f_transition[i - 1];

    f[i] = f[i - 1] + h * derivative_filter_mean[i - 1] + innovation;
    derivative_filter_mean[i] = derivative_filter_mean[i - 1]
      + cross_covariance / prediction_variance * innovation;
    derivative_filter_variance[i] = derivative_variance + lambda_i * h
      - square(cross_covariance) / prediction_variance;
  }
}
model {
  log_sigma_raw ~ std_normal();
  z_f_initial ~ std_normal();
  z_f_transition ~ std_normal();
  // SMOOTHING_PRIOR
  if (!prior_only)
    for (n in 1:N_obs)
      y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
}
generated quantities {
  vector[T] fprime;

  fprime[T] = normal_rng(
    derivative_filter_mean[T],
    sqrt(fmax(derivative_filter_variance[T], 0))
  );
  for (reverse_index in 1:(T - 1)) {
    int j = T - reverse_index;
    real h = deltat[j];
    real lambda_j = square(tau);
    real previous_variance = derivative_filter_variance[j];
    real prediction_variance =
      square(h) * previous_variance + lambda_j * h * square(h) / 3;
    real cross_next =
      h * previous_variance + lambda_j * square(h) / 2;
    real residual = f[j + 1] - f[j] - h * derivative_filter_mean[j];
    real previous_mean_given_f = derivative_filter_mean[j]
      + h * previous_variance / prediction_variance * residual;
    real previous_variance_given_f = previous_variance
      - square(h * previous_variance) / prediction_variance;
    real backward_covariance = previous_variance
      - h * previous_variance * cross_next / prediction_variance;
    real backward_gain =
      backward_covariance / derivative_filter_variance[j + 1];
    real backward_mean = previous_mean_given_f + backward_gain
      * (fprime[j + 1] - derivative_filter_mean[j + 1]);
    real backward_variance = previous_variance_given_f
      - square(backward_covariance) / derivative_filter_variance[j + 1];

    fprime[j] = normal_rng(
      backward_mean, sqrt(fmax(backward_variance, 0))
    );
  }
}
