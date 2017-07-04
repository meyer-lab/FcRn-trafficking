functions {
	real[] fcrn_model(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {
    real dC_dt[3];
    real Vc;
    real Vp;
    real Cl;
    real Q;
    real Qu;
    real sortF;
    real releaseF;
    real Vin;
    
    Vc = theta[1]; // Volume of the central compartment
    Vp = theta[2]; // Volume of the peripheral compartment
    Cl = theta[3];
    Q = theta[4];
    Qu = theta[5];
    sortF = theta[6];
    releaseF = theta[7];
    Vin = theta[8];

    // Cc, Cp, Cin
    dC_dt[1] = Q*(C[2] - C[1])/Vc - Cl*C[1]/Vc;
    dC_dt[2] = Q*(C[1] - C[2])/Vp - Qu*C[2]/Vp + Qu*C[3]*sortF*releaseF/Vp;
    dC_dt[3] = Qu*C[2]/Vin + Qu*C[3]/Vin*((1 - releaseF)*sortF - 1);

    return dC_dt;
	}

  real halfl_fcrn(real[] ts, real C0, real[] theta, real[] x_r, int[] x_i, real[] data_c, real varr) {
    real C_t0[3];
    real C_hat[6, 3];
    real sqErr;

    C_t0[1] = C0;
    C_t0[2] = 0.0;
    C_t0[3] = 0.0;

    C_hat = integrate_ode_rk45(fcrn_model, C_t0, 0.0, ts, theta, x_r, x_i);
  
    sqErr = dot_self(to_vector(C_hat[, 1]) - to_vector(data_c)) / varr; // Setup error measurement

    return sqErr;
  }
}
data {
  real wt_c[6];
  real ls_c[6];
  real dhs_c[6];
  real ts[6];
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
  real<lower=0> Vin; // Volume inside cells
  real<lower=0,upper=1> sortF_wt; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_wt; // Sorting fraction for release
  real<lower=0,upper=1> sortF_ls; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_ls; // Sorting fraction for release
  real<lower=0,upper=1> sortF_dhs; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_dhs; // Sorting fraction for release
  real<lower=0> varr; // Variance paramter
}
model {
  real theta[8];
  real sqErr;
  
  Vc ~ lognormal(0, 1.0);
  Vp ~ lognormal(0, 1.0);
  Cl ~ lognormal(-3, 1.0);
  Q ~ lognormal(0, 1.0);
  Qu ~ lognormal(-3, 1.0);
  varr ~ lognormal(0, 1.0);
  Vin ~ lognormal(0, 1.0);

  theta[1] = Vc;
  theta[2] = Vp;
  theta[3] = Cl;
  theta[4] = Q;
  theta[5] = Qu;
  theta[6] = sortF_wt;
  theta[7] = releaseF_wt;
  theta[8] = Vin;
  
  // Calculate data for wt condition
  sqErr = halfl_fcrn(ts, 14.0, theta, x_r, x_i, wt_c, varr);
  sqErr ~ chi_square(6); // Match to Student's t distribution
  
  
  // Calculate data for ls condition
  theta[6] = sortF_ls;
  theta[7] = releaseF_ls;
  
  sqErr = halfl_fcrn(ts, 24.0, theta, x_r, x_i, ls_c, varr);
  sqErr ~ chi_square(6); // Match to Student's t distribution
  
  
  // Calculate data for dhs condition
  theta[6] = sortF_dhs;
  theta[7] = releaseF_dhs;
  
  sqErr = halfl_fcrn(ts, 20.8, theta, x_r, x_i, dhs_c, varr);
  sqErr ~ chi_square(6); // Match to Student's t distribution
}
