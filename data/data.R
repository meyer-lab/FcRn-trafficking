# TODO: Spread data is stdev — check this is right
# TODO: Need values for FcRn KO

dataBase <- list(halflData = c(49.3, 335.9, 106.9, 204.3, 24),
	halflStd = c(2.7, 14.9, 4.3, 5.2, 1.0),
	ts = seq(1, 1000, length=100))

dataScarlett <- list(halflData = c(101.1, 323.0, 284.1, 294.7, 24),
	halflStd = c(11.4, 24.1, 9.7, 17.8, 1.0),
	ts = seq(1, 1000, length=100))


runsample <- function(name) {
	library(rstan);

	if (name == "diff") {
		dataIn <- dataBase
	} else if (name == "humanized") {
		dataIn <- dataScarlett
	}

	fit <- stan("model/diff.stan",
		cores = parallel::detectCores(),
		data = dataIn,
		iter = 500,
		chains = parallel::detectCores(),
		verbose = T,
		control = list(adapt_delta = 0.99));

	save(fit, file = paste("model/", name, "_samples.rds", sep = ""))
}
