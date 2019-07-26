
functions {
  real[] reduced_system (real time,
                         real[] y,
                         real[] yPK,
                         real[] theta,
                         real[] x_r, int[] x_i) { }
}

data {

}

transformed data {
  int nTheta = 8;
  int nCmt = 8;
  int nOde = 5;

  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);

  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_num_steps = 1000000;
  
  // . . .
}

parameters {

}

transformed parameters {
  // x = mixOde2CptModel_rk45(reduced_system, nOde,
  //                          time, amt, rate, ii, evid, cmt, addl, ss,
  //                          theta, biovar, tlag,
  //                          rel_tol, abs_tol, max_num_steps);
}

model {

}

generated quantities {

}
