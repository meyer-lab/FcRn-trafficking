functions {
  real[] fcrn_model(real t, real[] C, real[] th, real[] x_r, int[] x_i) {
    real dC_dt[3];
    real Cl;
    real Q;
    real Qu;
    real sortF;
    real releaseF;
    
    Q = th[2]; // Flow exchange rate between compartments
    Qu = th[3]; // Flow uptake into cells
    sortF = th[4]; // Sorting fraction for recycling
    releaseF = th[5]; // Sorting fraction for release

    // Cc, Cp, Cin
    dC_dt[1] = (Q*C[2] - Q*C[1])/1.0; // Vc is 1, no nonspecific clearance
    dC_dt[2] = (Q*C[1] - Q*C[2] - Qu*C[2] + Qu*C[3]*sortF*releaseF)/th[1]; // th[1] is Vp
    dC_dt[3] = Qu/th[6]*(C[2] + C[3]*((1 - releaseF)*sortF - 1));
    // th[6] is Vin

    return dC_dt;
  }
  real halfl_fcrn(real[] ts, real[] th, real[] x_r, int[] x_i) {
    real C_t0[3];
    real C_hat[1000, 3];
    real C0;

    C0 = 20.0; // TODO: What is the C0 concentration?

    C_t0[1] = C0;
    C_t0[2] = 0.0;
    C_t0[3] = 0.0;

    C_hat = integrate_ode_bdf(fcrn_model, C_t0, 0.0, ts, th, x_r, x_i);

    for (N in 1:size(ts)) {
      if (C_hat[N, 1] < C0/2) {
        return ts[N];
      }
    }

    return max(ts);
  }
}
data {
  real ts[1000];
}
transformed data {
  int x_i[1];
  real x_r[1];

  x_i[1] = 1;
  x_r[1] = 0.0;
}
parameters {
  real<lower=0> Vp; // The volume of the periferal compartment
  real<lower=0> Q; // Flow exchange rate between compartments
  real<lower=0> Qu; // Flow uptake into cells
  real<lower=0,upper=1> sortF_wt; // Sorting fraction for recycling
  real<lower=0,upper=1> sortF_ls; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_ls; // Sorting fraction for release
  real<lower=0,upper=1> sortF_yte; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_yte; // Sorting fraction for release
  real<lower=0,upper=1> sortF_dhs; // Sorting fraction for recycling
  real<lower=0> Vin; // Volume inside cells
}
transformed parameters {
  real actual_sortF_wt;
  real actual_sortF_dhs;
  real actual_sortF_ls;
  real actual_release_ls;
  actual_sortF_wt = sortF_wt * sortF_dhs * sortF_ls * sortF_yte;
  actual_sortF_dhs = sortF_dhs * sortF_ls * sortF_yte;
  actual_sortF_ls = sortF_ls * sortF_yte;
  actual_release_ls = releaseF_ls * releaseF_yte;
}
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
