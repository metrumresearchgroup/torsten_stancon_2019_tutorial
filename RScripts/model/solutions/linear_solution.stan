// Serves as a solution for the first exercise.

data {
  int N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real beta;
  real beta_0;
  real<lower = 0> sigma;
  // real log_sigma;
}

model {
  // priors
  beta ~ normal(2, 1);
  beta_0 ~ normal(25, 5);
  sigma ~ gamma(1, 1);
  
  // likelihood
  y ~ normal(beta * x + beta_0, sigma);
  // target += normal_lpdf(y | beta * x, sigma);
}

generated quantities {
  vector[N] y_pred = to_vector(normal_rng(beta * x + beta_0, sigma));
}

