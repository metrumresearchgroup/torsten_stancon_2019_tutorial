
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  
  // Event schedule
  int<lower = 1> cmt[nEvent];
  int evid[nEvent];
  int addl[nEvent];
  int ss[nEvent];
  real amt[nEvent];
  real time[nEvent];
  real rate[nEvent];
  real ii[nEvent];
  
  vector<lower = 0>[nObs] cObs;
}

transformed data {
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;
  int nCmt = 3;
  
  real biovar[nCmt] = {1.0, 1.0, 1.0};
  real tlag[nCmt] = {0.0, 0.0, 0.0};
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[nTheta] = {CL, Q, VC, VP, ka};
  row_vector<lower = 0>[nEvent] concentration;
  row_vector<lower = 0>[nObs] concentrationObs;
  matrix<lower = 0>[nCmt, nEvent] mass;

  mass = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss,
                          theta, biovar, tlag);
  // mass = PKModelTwoCpt(time, amt, rate, ii, evid, cmt, addl, ss,
  //                      theta, biovar, tlag);
  concentration = mass[2, ] ./ VC;
  concentrationObs = concentration[iObs];
}

model {
  // priors
  CL ~ lognormal(log(10), 0.25); 
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.25);
  VP ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(concentrationObs), sigma);
}

generated quantities {
  real concentrationObsPred[nObs];

  for (i in 1:nObs) {
    concentrationObsPred[i] = exp(normal_rng(log(concentrationObs[i]), sigma));
  }
}
