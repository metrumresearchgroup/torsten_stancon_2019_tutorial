functions{
    real[] oneCptPNODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    /* ... */
    return dxdt;
  }
}

data{
  int<lower = 1> nId;    /* population size */
  int<lower = 1> nt;            /* population's total time steps */
  int<lower = 1> nPNObs;
  int<lower = 1> nPNCens;
  int<lower = 1> iPNObs[nPNObs];
  int<lower = 1> iPNCens[nPNCens];
  real<lower = 0> time[nt];
  real<lower = 0> CL[nId];
  real<lower = 0> V[nId];

  /*
    create:
    pmx_solve_group_rk45() input args:
    amt, rate, ii, addl, cmt, evid, start, end, len, ss,
    nCmt, biovar(F), tlag
  */
}

parameters{
  real<lower = 0> ke0;
  real<lower = 0> alpha;
  real<lower = 0> beta;
}

transformed parameters{
  vector<lower = 0>[nPNObs] survObs;
  row_vector<lower = 0>[nPNObs] EdrugObs;
  vector<lower = 0>[nPNObs] hazardObs;
  vector<lower = 0>[nPNCens] survCens;
  matrix<lower = 0>[3, nt] x;

  /*
    create pmx_solve_group_rk45() argument "parms",
    use ODE function as hint.
  */
  
  /* 
     x = pmx_solve_group_rk45( ? ... );
  */

  for(i in 1:nPNObs) {
    survObs[i] = /* ? */    
  }
  EdrugObs = /* ? */ ;
  for(i in 1:nPNObs) {
    hazardObs[i] = /* ? */ ;
  }
  for(i in 1:nPNCens) {
    survCens[i] = /* ? */ ;
  }
}

model{
  ke0 ~ normal(0, 0.0005);
  alpha ~ normal(0, 0.000003);
  beta ~ normal(0, 1.5);
  
  target += log(hazardObs .* survObs); // observed PN event log likelihood
  target += log(survCens); // censored PN event log likelihood
}
