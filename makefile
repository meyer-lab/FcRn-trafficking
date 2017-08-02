
.PHONY: all shiny clean model sample

all: model/diff.rds

sample: model/diff_samples.rds model/humanized_samples.rds

clean:
	rm -f Analysis.pdf

model/diff.rds: model/diff.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

model/%_samples.rds: model/diff.rds
	R -e "source(\"data/data.R\"); runsample(\"$*\")"

Analysis.pdf: Analysis.Rmd
	Rscript -e "library(utils); rmarkdown::render('Analysis.Rmd', 'pdf_document')"

shiny: samples.rds
	R -e 'load("samples.rds"); shinystan::launch_shinystan(fit)'

predictions.rds: samples.rds
	R -e 'source("calcPredictions.R")'
