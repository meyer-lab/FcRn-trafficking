ddata <- list(ts = seq(1, 1000, length=1000)) # timepoints

runsample <- function() {
	library(rstan);

	fit <- stan("model/diff.stan", cores = parallel::detectCores(), data = ddata, chains = 12, verbose = T, control = list(adapt_delta = 0.99));

	save(fit, file = "samples.rds")
}
