
.PHONY: all shiny clean

all: model/diff.rds Analysis.pdf

clean:
	rm -f Analysis.pdf

model/%.rds: model/%.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

model/%_samples.rds: model/%.rds
	R -e "source(\"data/data.R\"); runsample(\"$*\")"

Analysis.pdf: Analysis.Rmd
	Rscript -e "library(utils); rmarkdown::render('Analysis.Rmd', 'pdf_document')"

shiny: samples.rds
	R -e 'load("samples.rds"); shinystan::launch_shinystan(fit)'

predictions.rds: samples.rds
	R -e 'source("calcPredictions.R")'
