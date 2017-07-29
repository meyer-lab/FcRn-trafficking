library(dplyr)
library(pbapply)
load('samples.rds')

expose_stan_functions("model/diff.stan")

samp_n <- 1000

sorts <- expand.grid(id = seq_len(samp_n),
                     sortF = seq(0, 0.98, length.out = 40),
                     releaseF = seq(0, 1, length.out = 3))

output <- as.data.frame(rstan::extract(fit)) %>%
  dplyr::select(Vp, Q, Qu, Vin) %>%
  dplyr::sample_n(samp_n) %>%
  dplyr::mutate(id = row_number()) %>%
  dplyr::full_join(sorts, by = c("id"))

# Make sure progress bar shows up
pboptions(type = "txt")

# Initiate cluster
cl <- parallel::makeCluster(parallel::detectCores() - 1)

# Run the predictions
output$halfl <- pbapply(output, 1, halfl_fcrn, cl = cl)

# Shutdown the cluster
parallel::stopCluster(cl)

save(output, file = 'predictions.rds', compress = "bzip2")
