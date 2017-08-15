
.PHONY: all test

all: 

test:
	R -e "library(methods); devtools::test('./fcrn')"

