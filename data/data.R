runsample <- function(name) {
	library(rstan);

	fit <- stan(paste("model/", name, ".stan", sep = ""),
		cores = parallel::detectCores(),
		data = list(ts = seq(1, 1000, length=1000)),
		chains = 12,
		verbose = T,
		control = list(adapt_delta = 0.99));

	save(fit, file = paste("model/", name, "_samples.rds", sep = ""))
}
