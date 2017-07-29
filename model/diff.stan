#include "common.stan"
model {
  real theta[6];
  real halfl;

  Vp ~ lognormal(0.0, 1.0);
  Q ~ lognormal(0.0, 1.0);
  Qu ~ lognormal(log(0.1), 0.5);
  Vin ~ lognormal(0, 1.0);

  theta[1] = Vp;
  theta[2] = Q;
  theta[3] = Qu;

  // wt recycles less than either engineered IgG's, so actual sorting is product
  theta[4] = actual_sortF_wt;
  theta[5] = 1.0;
  theta[6] = Vin;
  
  // Calculate data for wt condition
  halfl = halfl_fcrn(ts, theta, x_r, x_i);

  halfl ~ normal(49.3, 2.7); // TODO: This is the stderr, check if correct


  // dhs recycles less than ls, so actual sorting is product
  theta[4] = actual_sortF_dhs;
  
  halfl = halfl_fcrn(ts, theta, x_r, x_i);

  halfl ~ normal(335.9, 14.9); // TODO: This is the stderr, check if correct
  

  // Calculate data for ls condition
  theta[4] = actual_sortF_ls;
  theta[5] = actual_release_ls;
  
  halfl = halfl_fcrn(ts, theta, x_r, x_i);

  halfl ~ normal(106.9, 4.3); // TODO: This is the stderr, check if correct


  // Calculate data for yte condition
  theta[4] = sortF_yte;
  theta[5] = releaseF_yte;
  
  halfl = halfl_fcrn(ts, theta, x_r, x_i);

  halfl ~ normal(204.3, 5.2); // TODO: This is the stderr, check if correct


  // Calculate data for FcRn KO
  theta[4] = 0.0;

  halfl = halfl_fcrn(ts, theta, x_r, x_i);

  halfl ~ normal(24, 1.0); // TODO: Need values for this
}
