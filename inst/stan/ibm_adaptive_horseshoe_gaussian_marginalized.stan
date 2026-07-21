// Gaussian-data, marginalized version of ibm_adaptive_horseshoe.stan.
//
// The two-dimensional IBM state is integrated out with a Kalman filter.  The
// parameters and priors below deliberately match the existing noncentered
// horseshoe model so that differences in sampling performance are attributable
// to marginalization rather than to a changed statistical model.

data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  vector[N_obs] y_obs;

  real log_sigma_mu;
  real<lower=0> log_sigma_sd;
  real<lower=0> zeta;
  real<lower=0> initial_sd;
}

parameters {
  real log_sigma;
  vector<lower=0, upper=1>[T - 1] z_tau;
  real<lower=0, upper=1> z_gamma;
}

transformed parameters {
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);
  real<lower=0> gamma = zeta * tan(z_gamma * pi() / 2);
  vector<lower=0>[T - 1] tau;

  for (i in 1:(T - 1))
    tau[i] = gamma * tan(z_tau[i] * pi() / 2);
}

model {
  // Filtering state is ordered (fprime, f).  On this ordering the correlated
  // unit-time IBM covariance is initial_sd^2 * [[1, 1/2], [1/2, 1/3]].
  vector[2] m = rep_vector(0, 2);
  matrix[2, 2] P = square(initial_sd) *
    [[1.0, 0.5], [0.5, 1.0 / 3.0]];

  log_sigma ~ std_normal();
  z_tau ~ uniform(0, 1);
  z_gamma ~ uniform(0, 1);

  for (i in 1:T) {
    if (i > 1) {
      real h = deltat[i - 1];
      matrix[2, 2] F = [[1.0, 0.0], [h, 1.0]];
      matrix[2, 2] S = [[h, 0.5 * square(h)],
                        [0.5 * square(h), h * square(h) / 3.0]];

      m = F * m;
      P = F * P * F' + square(tau[i - 1]) * S;
      P = 0.5 * (P + P');
    }

    // Sequential scalar updates handle replicate observations at a grid point.
    // obs_time_idx is nondecreasing in data prepared by ibm_smooth(); the
    // explicit scan also makes the model robust to infer_at-only grid points.
    for (n in 1:N_obs) {
      if (obs_time_idx[n] == i) {
        real innovation = y_obs[n] - m[2];
        real innovation_var = P[2, 2] + square(sigma);
        vector[2] gain = P[, 2] / innovation_var;
        row_vector[2] level_cov = P[2, ];

        target += normal_lpdf(innovation | 0, sqrt(innovation_var));
        m += gain * innovation;
        P -= gain * level_cov;
        P = 0.5 * (P + P');
      }
    }
  }
}
