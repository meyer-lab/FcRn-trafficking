library(dplyr)
load('samples.rds')

halfl_fcrn <- function(th) {
  fcrn_model <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dCc = Q*Cp - Q*Cc
      dCp = (Q*Cc - Q*Cp - Qu*Cp + Qu*Cin*sortF*releaseF)/Vp
      dCin = Qu/Vin*(Cp + Cin*((1 - releaseF)*sortF - 1))
      
      list(c(dCc, dCp, dCin))
    })
  }
  
  C0 <- 20.0
  
  if (th['sortF']*th['releaseF'] < 0.4) {
    ts <- seq_len(100)
  } else {
    ts <- seq(0, 10000, length = 3000)
  }
  
  odeSol <- deSolve::ode(y = c(Cc = C0, Cp = 0, Cin = 0),
                         times = ts, func = fcrn_model, parms = th)
  
  return(approx(odeSol[,2], odeSol[,1], C0/2.0)$y)
}

samp_n <- 1000

sorts <- expand.grid(id = seq_len(samp_n),
                     sortF = seq(0, 0.95, length.out = 20),
                     releaseF = seq(0, 1, length.out = 3))

output <- as.data.frame(rstan::extract(fit)) %>%
  dplyr::select(Vp, Q, Qu, Vin) %>%
  dplyr::sample_n(samp_n) %>%
  dplyr::mutate(id = row_number()) %>%
  dplyr::full_join(sorts, by = c("id"))

# Initiate cluster
cl <- parallel::makeCluster(parallel::detectCores() - 1)

# Run the predictions
output$halfl <- pbapply::pbapply(output, 1, halfl_fcrn, cl = cl)

# Shutdown the cluster
parallel::stopCluster(cl)

save(output, file = 'predictions.rds')
