
.PHONY: all shiny clean model sample

all: model/diff.rds Analysis.pdf ModelDescription.pdf

sample: model/diff_samples.rds model/humanized_samples.rds

clean:
	rm -f Analysis.pdf

model/diff.rds: model/diff.stan
	R -e 'rstan::stan_model("model/diff.stan", auto_write = T)'

model/%_samples.rds: model/diff.rds
	R -e "source(\"data/data.R\"); runsample(\"$*\")"

%.pdf: %.Rmd
	Rscript -e "library(utils); rmarkdown::render('$*.Rmd', 'pdf_document')"

diff_shiny: model/diff_samples.rds
	R -e 'load("model/diff_samples.rds"); shinystan::launch_shinystan(fit)'

humanized_shiny: model/humanized_samples.rds
	R -e 'load("model/humanized_samples.rds"); shinystan::launch_shinystan(fit)'

predictions.rds: samples.rds
	R -e 'source("calcPredictions.R")'
