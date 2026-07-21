functions {
  real log_gn_density(real f, real fprime, 
                      real f_prev, real fprime_prev,
                      real gamma, real alpha, 
                      matrix L_t,
                      real dt) {
    vector[2] diff;

    // State ordering: [f', f]
    diff[1] = fprime - fprime_prev;
    diff[2] = f      - (f_prev + fprime_prev * dt);

    // z = L_t^{-1} * diff => z'z = diff' * K^{-1} * diff
    vector[2] z = mdivide_left_tri_low(L_t, diff);
    real quadform = dot_self(z);  // = diff' K^{-1} diff

    // dimension d = 2 => constant contributes -2 * log(gamma)
    // generalized normal radial part has (quadform / gamma^2)^(alpha/2)
    real t1 = -2 * log(gamma);
    real t2 = -0.5 * pow(quadform / square(gamma), alpha / 2);

    return t1 + t2;  // plus constants independent of parameters (dropped)
  }
}

data {
  int<lower=1> N_obs; // number of observations (including replicates if there are any)
  int<lower=1> T; // number of latent grid points (unique times + infer_at)
  array[N_obs] int<lower=1, upper=T> obs_time_idx;  // maps each obs to a grid index
  vector[T-1] deltat; // time diffs on the (scaled) grid
  vector[N_obs] y_obs; // observed y values (also scaled)
  int<lower=0, upper=1> regular;

  // gaussian likelihood - lognormal prior for sigma
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;

  // ibm - global scale for tau/gamma
  real zeta;
  // bridge exponent
  real<lower=0, upper=2> alpha;
}

transformed data{
  array[T-1] matrix[2,2] Lt;

  Lt[1] = to_matrix( // cholesky decomp of ibm covariance
    [[sqrt(deltat[1]), 0],
    [(deltat[1]^1.5) * 0.5, (deltat[1]^1.5) * inv_sqrt(12)]]);

  for(i in 2:(T-1)) {
    if (regular == 0) {
      Lt[i] = to_matrix( // cholesky decomp of ibm covariance
        [[sqrt(deltat[i]), 0],
        [(deltat[i]^1.5) * 0.5, (deltat[i]^1.5) * inv_sqrt(12)]]);
    } else {
      Lt[i] = Lt[1];
    }
  }
}

parameters {
  real log_sigma;
  vector[T] f;
  vector[T] fprime;
  real<lower=0, upper=1> z_gamma;

}

transformed parameters {
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);
  real gamma = zeta * tan(z_gamma * pi() / 2); // half-Cauchy(0, zeta)
}

model {
  // global shrinkage parameter
  z_gamma ~ uniform(0, 1);

  // hyperpriors
  log_sigma ~ std_normal();


  // ibm prior
  fprime[1] ~ std_normal();
  f[1] ~ std_normal();
  for (i in 2:T) {
    target += log_gn_density(f[i], fprime[i],
                             f[i-1], fprime[i-1],
                             gamma, alpha, Lt[i-1],
                             deltat[i-1]);
  }

  // likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
