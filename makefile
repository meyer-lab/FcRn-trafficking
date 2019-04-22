.PHONY: all clean

all: output.pdf

output.pdf:
	R CMD INSTALL --no-multiarch --with-keep.source fcrn
	R -e "library(fcrn); fcrn::full_plot()"

clean:
	rm -f output.pdf Rplots.pdf
