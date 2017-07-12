
.PHONY: all clean shiny

all: model/diff.rds

clean:
	rm -f model/diff.rds samples.rds

model/%.rds: model/%.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

samples.rds: model/diff.stan model/diff.rds
	R -e 'library(rstan); options(mc.cores = parallel::detectCores()); source("data/data.R"); fit <- stan("model/diff.stan", data = ddata, chains = 8, verbose = T, iter = 500, control = list(adapt_delta = 0.99)); save(fit, file = "samples.rds")'

shiny: samples.rds
	R -e 'load("samples.rds"); shinystan::launch_shinystan(fit)'
