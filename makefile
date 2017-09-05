
.PHONY: all clean

all: output.pdf README.pdf README.docx

output.pdf:
	R CMD INSTALL --no-multiarch --with-keep.source fcrn
	R -e "library(fcrn); fcrn::full_plot()"

README.pdf: README.md
	pandoc -o README.pdf README.md

README.docx: README.md
	pandoc -o README.docx README.md

clean:
	rm -f output.pdf README.pdf README.docx Rplots.pdf