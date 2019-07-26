
functions {
  real[] system (real time,
                 real[] y,
                 real[] theta,
                 real[] x_r, int[] x_i) {
    real CL = theta[1];
    real Q = theta[2];
    real VC = theta[3];
    real VP = theta[4];
    real ka = theta[5];
    real mtt = theta[6];
    real ktr = 4 / mtt;
    real circ0 = theta[7];
    real gamma = theta[8];

    real alpha = 3e-4;
    real Edrug;

    // Device to implement initial condition
    real prol = y[4] + circ0;
    real transit1 = y[5] + circ0;
    real transit2 = y[6] + circ0;
    real transit3 = y[7] + circ0;
    real circ = fmax(machine_precision(), y[8] + circ0);
    // real circ = y[8] + circ0;

    real dydt[8];

    dydt[1] = -ka * y[1];
    dydt[2] = ka * y[1] - (CL + Q) / VC * y[2] + Q / VP * y[3];
    dydt[3] = Q / VC * y[2] - Q / VP * y[3];

    Edrug = fmin(1.0, alpha * y[2] / VC);
    // Edrug = alpha * y[2] / VC;

    dydt[4] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
    dydt[5] = ktr * (prol - transit1);
    dydt[6] = ktr * (transit1 - transit2);
    dydt[7] = ktr * (transit2 - transit3);
    dydt[8] = ktr * (transit3 - circ);

    return dydt;
  }
}

data {
  int<lower = 1> nEvent;
  int<lower = 1> nObsPK;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPD[nObsPD];

  // Event schedule
  int<lower = 1> cmt[nEvent];
  int evid[nEvent];
  int addl[nEvent];
  int ss[nEvent];
  real amt[nEvent];
  real time[nEvent];
  real rate[nEvent];
  real ii[nEvent];

  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
}

transformed data {
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nTheta = 8;
  int nCmt = 8;

  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);

  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  // int max_num_steps = 1000000;

  // for testing purposes
  // real CL_x = 10;
  // real Q_x = 15;
  // real VC_x = 35;
  // real VP_x = 105;
  // real ka_x = 2;
  // real mtt_x = 125;
  // real circ0_x = 5.0;
  // real gamma_x = 0.17;
  // real alpha_x = 3e-4;
  // real theta_x[nTheta] = {CL_x, Q_x, VC_x, VP_x, ka_x,
  //                         mtt_x, circ0_x, gamma_x};
}

parameters {
  // PK parameters
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;

  // PD parameters
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> gamma;
  real<lower = 0> sigmaNeut;
}

transformed parameters {
  real theta[nTheta] = {CL, Q, VC, VP, ka,
                        mtt, circ0, gamma};

  vector<lower = 0>[nObsPK] concentrationObs;
  vector<lower = 0>[nObsPD] neut;

  // matrix to store drug mass and neutrophil count.
  matrix[nEvent, nCmt] x;

  x = generalOdeModel_rk45(system, nCmt,
                          time, amt, rate, ii, evid, cmt, addl, ss,
                          theta, biovar, tlag,
                          rel_tol, abs_tol, 1e+5);

  concentrationObs = x[iObsPK, 2] ./ VC;

  neut = x[iObsPD, 8] + circ0;
}

model {
  // priors for PK parameters
  CL ~ lognormal(log(10), 0.25); 
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.25);
  VP ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);

  // priors for PD parameters
  mtt ~ lognormal(log(125), 0.2);
  circ0 ~ lognormal(log(5), 0.2);
  gamma ~ lognormal(log(0.17), 0.1);
  sigmaNeut ~ cauchy(0, 1);

  logCObs ~ normal(log(concentrationObs), sigma);
  logNeutObs ~ normal(log(neut), sigmaNeut);
}

generated quantities {
  real concentrationObsPred[nObsPK];
  real neutObsPred[nObsPD];

  for (i in 1:nObsPK)
    concentrationObsPred[i] = exp(normal_rng(log(concentrationObs[i]), sigma));

  for (i in 1:nObsPD)
    neutObsPred[i] = exp(normal_rng(log(neut[i]), sigmaNeut));
}
