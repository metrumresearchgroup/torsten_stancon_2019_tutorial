
functions {
  real[] system(real time, real[] y, 
                real[] theta, real[] x_r, int[] x_i) {
  real dydt[3];
  real CL = theta[1];
  real Q = theta[2];
  real VC = theta[3];
  real VP = theta[4];
  real ka = theta[5];
  
  dydt[1] = - ka * y[1];
  dydt[2] = ka * y[1] - (CL + Q) / VC * y[2] + Q / VP * y[3];
  dydt[3] = Q / VC * y[2] - Q / VP * y[3];

  return dydt;
  }
}

data {
  int nEvent;
  int nObs;
  int iObs[nObs];
  
  // Event schedule
  real time[nEvent];
  real amt[nEvent];
  real rate[nEvent];
  real ii[nEvent];
  int ss[nEvent];
  int evid[nEvent];
  int cmt[nEvent];
  int addl[nEvent];

  // observation
  vector[nObs] cObs;
}

transformed data {
  int nCmt = 3;
  int nParm = 5;
  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);
  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_num_steps = 100000;
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real <lower = 0> sigma;
}

transformed parameters {
  vector[nObs] concentration;
  matrix[nEvent, nCmt] mass;
  real theta[nParm] = {CL, Q, VC, VP, ka};

  mass = generalOdeModel_rk45(system, nCmt,
                              time, amt, rate, ii, evid, cmt,
                              addl, ss,
                              theta, biovar, tlag,
                              rel_tol, abs_tol, max_num_steps);

  // mass = PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
  //                      theta, biovar, tlag);

  concentration = mass[iObs, 2] / VC;
}

model {
  // prior
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.25);
  VP ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);

  // likelihood
  cObs ~ lognormal(log(concentration), sigma);
}

generated quantities {
  vector[nObs] concentrationObsPred;

  for (i in 1:nObs)
    concentrationObsPred[i] 
      = lognormal_rng(log(concentration[i]), sigma);
}
