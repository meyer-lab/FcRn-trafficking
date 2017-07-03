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

  vector halfl_fcrn(real[] ts, real C0, real[] theta, real[] x_r, int[] x_i) {
    real C_t0[2];
    real C_hat[6, 2];

    C_t0[1] = C0;
    C_t0[2] = 0.0;

    C_hat = integrate_ode_rk45(fcrn_model, C_t0, 0.0, ts, theta, x_r, x_i);

    return to_vector(C_hat[, 1]);
  }
}
data {
  real<lower=0> ts[6];
  real wt_c[6];
}
transformed data {
  real x_r[1];
  int x_i[1];
  
  x_r[1] = 0.0;
  x_i[1] = 1;
}
parameters {
  real<lower=0> Vc; // The volume of the central compartment
  real<lower=0> Vp; // The volume of the periferal compartment
  real<lower=0> Cl; // Nonspecific clearance rate
  real<lower=0> Q; // Flow exchange rate between compartments
  real<lower=0> Qu; // Flow uptake into cells
  real<lower=0> varr; // Variance paramter
}
model {
  real theta[5];
  vector[6] c_hat;
  real sqErr;
  
  Vc ~ lognormal(0, 1.0);
  Vp ~ lognormal(0, 1.0);
  Cl ~ lognormal(0, 1.0);
  Q ~ lognormal(0, 1.0);
  Qu ~ lognormal(0, 1.0);
  varr ~ lognormal(0, 1.0);

  theta[1] = Vc;
  theta[2] = Vp;
  theta[3] = Cl;
  theta[4] = Q;
  theta[5] = Qu;
  
  c_hat = halfl_fcrn(ts, 14.0, theta, x_r, x_i);
  
  c_hat[6] = wt_c[6];
  
  //sqErr = dot_self(log(c_hat) - log(to_vector(wt_c)));
  
  //sqErr ~ student_t(5, 0.0, varr);
  
  
}
