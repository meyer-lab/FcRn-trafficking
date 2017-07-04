functions {
	real[] fcrn_model(real t, real[] C, real[] th, real[] x_r, int[] x_i) {
    real dC_dt[3];
    real Cl;
    real Q;
    real Qu;
    real sortF;
    real releaseF;
    
    Cl = th[3]; // Nonspecific clearance rate
    Q = th[4]; // Flow exchange rate between compartments
    Qu = th[5]; // Flow uptake into cells
    sortF = th[6]; // Sorting fraction for recycling
    releaseF = th[7]; // Sorting fraction for release

    // Cc, Cp, Cin
    dC_dt[1] = (Q*C[2] - Q*C[1] - Cl*C[1])/th[1]; // th[1] is Vc
    dC_dt[2] = (Q*C[1] - Q*C[2] - Qu*C[2] + Qu*C[3]*sortF*releaseF)/th[2]; // th[2] is Vp
    dC_dt[3] = Qu/th[8]*(C[2] + C[3]*((1 - releaseF)*sortF - 1));
    // th[8] is Vin

    return dC_dt;
	}

  real halfl_fcrn(real[] ts, real C0, real[] th, real[] x_r, int[] x_i, real[] data_c) {
    real C_t0[3];
    real C_hat[6, 3];

    C_t0[1] = C0;
    C_t0[2] = 0.0;
    C_t0[3] = 0.0;

    C_hat = integrate_ode_bdf(fcrn_model, C_t0, 0.0, ts, th, x_r, x_i);

    return dot_self(to_vector(C_hat[, 1]) ./ to_vector(data_c)); // Setup error measurement
  }
}
data {
  real wt_c[6];
  real ls_c[6];
  real dhs_c[6];
  real ko_c[6];
  real ts[6];
}
transformed data {
  real x_r[1];
  int x_i[1];
  
  x_r[1] = 0.0;
  x_i[1] = 1;
}
parameters {
  real<lower=0,upper=100> Vc; // The volume of the central compartment
  real<lower=0,upper=100> Vp; // The volume of the periferal compartment
  real<lower=0,upper=1000> Q; // Flow exchange rate between compartments
  real<lower=0,upper=1000> Qu; // Flow uptake into cells
  real<lower=0,upper=1> sortF_wt; // Sorting fraction for recycling
  real<lower=0,upper=1> sortF_ls; // Sorting fraction for recycling
  real<lower=0,upper=1> releaseF_ls; // Sorting fraction for release
  real<lower=0,upper=1> sortF_dhs; // Sorting fraction for recycling
  real<lower=0,upper=100> Vin; // Volume inside cells
  real<lower=0,upper=10> varr; // Variance paramter
}
model {
  real theta[8];
  real sqErr;
  
  Vc ~ lognormal(0, 0.5);
  Vp ~ lognormal(0, 1.0);
  Q ~ lognormal(0, 1.0);
  Qu ~ lognormal(-2, 1.0);
  varr ~ lognormal(-2, 0.5);
  Vin ~ lognormal(0, 1.0);

  theta[1] = Vc;
  theta[2] = Vp;
  theta[3] = 0.0; // Set nonspecific clearance to 0
  theta[4] = Q;
  theta[5] = Qu;
  theta[6] = sortF_wt;
  theta[7] = 1.0;
  theta[8] = Vin;
  
  // Calculate data for wt condition
  sqErr = halfl_fcrn(ts, 14.0, theta, x_r, x_i, wt_c) / varr;
  sqErr ~ chi_square(6); // Match to Student's t distribution


  // Calculate data for dhs condition
  theta[6] = sortF_dhs;
  
  sqErr = halfl_fcrn(ts, 20.8, theta, x_r, x_i, dhs_c) / varr;
  sqErr ~ chi_square(6); // Match to Student's t distribution
  
  
  // Calculate data for ls condition
  theta[6] = sortF_ls;
  theta[7] = releaseF_ls;
  
  sqErr = halfl_fcrn(ts, 24.0, theta, x_r, x_i, ls_c) / varr;
  sqErr ~ chi_square(6); // Match to Student's t distribution


  // Calculate data for FcRn KO
  theta[6] = 0.0;

  sqErr = halfl_fcrn(ts, 20.0, theta, x_r, x_i, ko_c) / varr;
  sqErr ~ chi_square(6); // Match to Student's t distribution
}
