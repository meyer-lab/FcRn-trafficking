functions {
  // Precise method, which guarantees v = v1 when t = 1.
  real lerp(real v0, real v1, real t) {
    return (1 - t) * v0 + t * v1;
  }
  real fcrn_model(real t, matrix A, vector C_t0) {
    return (matrix_exp(t * A) * C_t0)[1];
  }
  matrix make_Matrix(real[] th) {
    real Cl;
    real Q;
    real Qu;
    real sortF;
    real releaseF;
    real Vp;
    real Vin;
    matrix[3, 3] A;
    
    Q = th[2]; // Flow exchange rate between compartments
    Qu = th[3]; // Flow uptake into cells
    sortF = th[4]; // Sorting fraction for recycling
    releaseF = th[5]; // Sorting fraction for release
    Vp = th[1];
    Vin = th[6];
    
    A[1, 1] = -Q;
    A[1, 2] = Q;
    A[1, 3] = 0;
    A[2, 1] = Q/Vp;
    A[2, 2] = -(Q+Qu)/Vp;
    A[2, 3] = Qu*sortF*releaseF/Vp;
    A[3, 1] = 0;
    A[3, 2] = Qu/Vin;
    A[3, 3] = Qu*((1-releaseF)*sortF-1)/Vin;
    
    return(A);
  }
  real halfl_fcrn(real[] ts, real[] th) {
    real C_hat[100];
    real inter;
    real C0;
    vector[3] C_t0;
    matrix[3, 3] A;
    
    A = make_Matrix(th);
    
    C0 = 20.0; // TODO: What is the C0 concentration?
    
    C_t0[1] = C0;
    C_t0[2] = 0.0;
    C_t0[3] = 0.0;
    
    for (N in 1:size(ts)) {
      C_hat[N] = fcrn_model(ts[N], A, C_t0);
    }

    for (N in 2:size(ts)) {
      if (C_hat[N] < C0/2) {
        // Find interpolation point
        inter = (C0/2 - C_hat[N-1]) / (C_hat[N] - C_hat[N-1]);

        return lerp(ts[N-1], ts[N], inter);
      }
    }

    return max(ts);
  }
}
data {
  real ts[100];
  real halflData[5];
  real halflStd[5];
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
  halfl = halfl_fcrn(ts, theta);

  halfl ~ normal(halflData[1], halflStd[1]);


  // dhs recycles less than ls, so actual sorting is product
  theta[4] = actual_sortF_dhs;
  
  halfl = halfl_fcrn(ts, theta);

  halfl ~ normal(halflData[2], halflStd[2]);
  

  // Calculate data for ls condition
  theta[4] = actual_sortF_ls;
  theta[5] = actual_release_ls;
  
  halfl = halfl_fcrn(ts, theta);

  halfl ~ normal(halflData[3], halflStd[3]);


  // Calculate data for yte condition
  theta[4] = sortF_yte;
  theta[5] = releaseF_yte;
  
  halfl = halfl_fcrn(ts, theta);

  halfl ~ normal(halflData[4], halflStd[4]);


  // Calculate data for FcRn KO
  theta[4] = 0.0;

  halfl = halfl_fcrn(ts, theta);

  halfl ~ normal(halflData[5], halflStd[5]);
}
