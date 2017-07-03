functions {
	real[] fcrn_model(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {
    real dC_dt[2];
    real Vc;
    real Vp;
    real Cl;
    real Q;
    real Qu;
    
    
    Vc = theta[1]; // Volume of the central compartment
    Vp = theta[2]; // Volume of the peripheral compartment
    Cl = theta[3];
    Q = theta[4];
    Qu = theta[5];

    // Cc, Cp
    dC_dt[1] = Q*(C[2] - C[1])/Vc - Cl*C[1]/Vc;
    dC_dt[2] = Q*(C[1] - C[2])/Vp - Qu*C[2]/Vp;

    return dC_dt;
	}

  real halfl_fcrn(int N_t, real[] ts, real C0, real[] theta, real[] x_r, int[] x_i) {
    real C_t0[2];
    real C_hat[N_t,2];

    C_t0[1] = C0;
    C_t0[2] = 0.0;

    C_hat = integrate_ode_rk45(fcrn_model, C_t0, 0.0, ts, theta, x_r, x_i);

    for (ii in 1:N_t) {
      if (C_hat[ii,1] < C0/2.0) {
        return ts[ii];
      }
    }

    return ts[N_t];
  }
}
data {
  int<lower=0> N_obs;
  real<lower=0> obs_t[N_obs]; // timepoint ids
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  real<lower=0> Vc; // The volume of the central compartment
  real<lower=0> Vp; // The volume of the periferal compartment
  real<lower=0> Cl; // Nonspecific clearance rate
  real<lower=0> Q; // Flow exchange rate between compartments
  real<lower=0> Qu; // Flow uptake into cells
}
model {
  real theta[5];
  real halfl;
  
  Vc ~ lognormal(0, 0.1);
  Vp ~ lognormal(0, 0.1);
  Cl ~ lognormal(0, 0.1);
  Q ~ lognormal(0, 0.1);
  Qu ~ lognormal(0, 0.1);

  theta[1] = Vc;
  theta[2] = Vp;
  theta[3] = Cl;
  theta[4] = Q;
  theta[5] = Qu;
  
  halfl = halfl_fcrn(N_obs, obs_t, 5.0, theta, x_r, x_i);
  halfl ~ normal(41.2, 10.0);
}
