
.PHONY: all shiny

all: model/diff.rds

model/%.rds: model/%.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

samples.rds: model/diff.stan model/diff.rds
	R -e 'source("data/data.R"); runsample()'

shiny: samples.rds
	R -e 'load("samples.rds"); shinystan::launch_shinystan(fit)'
