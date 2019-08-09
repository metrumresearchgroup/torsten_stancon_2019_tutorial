functions{
  real[] reaction(real t, real[] x, real[] theta, real[] r, int[] i){
    real dxdt[3];
    real k1 = theta[1];
    real k2 = theta[2];
    real k3 = theta[3];
    dxdt[1] = -k1*x[1] + k3*x[2]*x[3];
    dxdt[2] =  k1*x[1] - k3*x[2]*x[3] - k2*(x[2])^2;
    dxdt[3] =  k2*(x[2])^2;
    return dxdt;
  }
}

data {
  int<lower=1> nsub;
  int<lower=1> len[nsub];  
  int<lower=1> ntot;  
  real ts[ntot];
  real obs[ntot];
}

transformed data {
  int i1[nsub];
  int i2[nsub];
  real t0 = 0.0;
  real xr[0];
  int xi[0];
  real theta[3] = {0.04, 3.0e7, 1.0e4};
  i1[1] = 1;
  i2[1] = len[1];
  for (i in 2:nsub) {
    i1[i] = i2[i-1] + 1;
    i2[i] = i1[i] + len[i] - 1;
  }
}

parameters {
  /*  p1=0.04, p2=1e4, and p3=3e7 */
  real<lower = 0> y0_mu;
  real<lower = 0> y0_1[nsub];
  real<lower = 0> sigma;
}

transformed parameters {
  real y0[3];
  real x[ntot, 3];
  real x3[ntot];
  for (i in 1:nsub) {
    y0[1] = y0_1[i];
    y0[2] = 0.0;
    y0[3] = 0.0;
    x[i1[i]:i2[i], ] = pmx_integrate_ode_bdf(reaction, y0, t0, ts[i1[i]:i2[i]], theta, xr, xi, 1.e-4, 1.e-8, 10000);
  }
  x3 = x[ , 3];
}

model {
  y0_mu ~ lognormal(log(2.0), 0.5);
  for (i in 1:nsub) {
    y0_1[i] ~ lognormal(y0_mu, 0.5);    
  }
  sigma ~ cauchy(0, 0.5); 
  obs ~ lognormal(log(x3), sigma);
}
