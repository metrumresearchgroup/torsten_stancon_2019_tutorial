
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  int<lower = 1> nSubjects;
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];

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
  int nIIV = 5;

  real biovar[nCmt] = {1.0, 1.0, 1.0};
  real tlag[nCmt] = {0.0, 0.0, 0.0};
}

parameters {
  // Population parameters
  real<lower = 0> CL_pop;
  real<lower = 0> Q_pop;
  real<lower = 0> VC_pop;
  real<lower = 0> VP_pop;
  real<lower = 0> ka_pop;

  // Inter-individual variability
  vector<lower = 0, upper = 2>[nIIV] omega;
  real<lower = 0> theta[nSubjects, nTheta];

  real<lower = 0> sigma;
}

transformed parameters {
  vector<lower = 0>[nTheta] 
    theta_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});
  vector<lower = 0>[nEvent] concentration;
  vector<lower = 0>[nObs] concentrationObs;
  matrix<lower = 0>[nEvent, nCmt] mass;

  for (j in 1:nSubjects) {
    mass[start[j]:end[j]] = PKModelTwoCpt(time[start[j]:end[j]],
                                          amt[start[j]:end[j]],
                                          rate[start[j]:end[j]],
                                          ii[start[j]:end[j]],
                                          evid[start[j]:end[j]],
                                          cmt[start[j]:end[j]],
                                          addl[start[j]:end[j]],
                                          ss[start[j]:end[j]],
                                          theta[j, ], biovar, tlag);

    concentration[start[j]:end[j]] = 
                      mass[start[j]:end[j], 2] / theta[j, 3];
  }

  concentrationObs = concentration[iObs];
}

model {
  // priors
  CL_pop ~ lognormal(log(10), 0.25); 
  Q_pop ~ lognormal(log(15), 0.5);
  VC_pop ~ lognormal(log(35), 0.25);
  VP_pop ~ lognormal(log(105), 0.5);
  ka_pop ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);
  omega ~ cauchy(0, 0.2);

  // hierarchical likelihood
  for (j in 1:nSubjects)
    theta[j, ] ~ lognormal(log(theta_pop), omega);

  logCObs ~ normal(log(concentrationObs), sigma);
}


generated quantities {
    real log_theta[nSubjects, nTheta] = log(theta);
}
