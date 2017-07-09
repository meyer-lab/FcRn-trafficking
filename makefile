
.PHONY: all clean

all: model/diff.rds

clean:
	rm -f model/diff.rds samples.rds

model/%.rds: model/%.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

samples.rds: model/diff.stan model/diff.rds
	R -e 'library(rstan); options(mc.cores = parallel::detectCores()); source("data/data.R"); fit <- stan("model/diff.stan", data = ddata, verbose = T); save(fit, file = "samples.rds")'