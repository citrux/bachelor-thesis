all: bachelor-thesis.pdf

%.pdf: %.tex
	latexmk -pdf $<
