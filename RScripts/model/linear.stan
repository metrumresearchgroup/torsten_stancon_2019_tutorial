// template for Stan file.

data {
  int<lower = 1> N;  // number of data points
  vector[N] x;
  vector[N] y;
}

parameters {
  real beta;
  // real beta_0;
  real<lower = 0> sigma;
}

model {
  // prior distribution
  beta ~ normal(2.0, 1.0);
  sigma ~ gamma(1.0, 1.0);

  // likelihood distribution
  y ~ normal(beta * x, sigma);
  // y ~ normal(beta_0 + beta * x, sigma);
}

generated quantities {
  vector[N] y_pred;
  
  for (i in 1:N)
    y_pred[i] = normal_rng(beta * x[i], sigma);
    // y_pred[i] = normal_rng(beta_0 + beta * x[i], sigma);
}
