
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
  int nParms = 5;
  int nCmt = 3;
  real biovar[3] = rep_array(1.0, nCmt);
  real tlag[3] = rep_array(0.0, nCmt);
}

parameters {
  real CL_pop;
  real Q_pop;
  real VC_pop;
  real VP_pop;
  real ka_pop;
  
  vector[nParms] omega;
  matrix<lower = 0>[nSubjects, nParms] theta;
  
  real<lower = 0> sigma;
}

transformed parameters {
  vector[nParms] theta_pop 
    = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});
  matrix[nEvent, nCmt] mass;
  vector[nEvent] concentration;
  vector[nObs] concentrationObs;

  for (j in 1:nSubjects) {
    mass[start[j]:end[j], ] = PKModelTwoCpt(time[start[j]:end[j]],
                                          amt[start[j]:end[j]],
                                          rate[start[j]:end[j]],
                                          ii[start[j]:end[j]],
                                          evid[start[j]:end[j]],
                                          cmt[start[j]:end[j]],
                                          addl[start[j]:end[j]],
                                          ss[start[j]:end[j]],
                                          to_array_1d(theta[j, ]), 
                                          biovar, tlag);

   concentration[start[j]:end[j]] 
     = mass[start[j]:end[j], 2] / theta[j, 3];
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

  cObs ~ lognormal(log(concentrationObs), sigma);
}
