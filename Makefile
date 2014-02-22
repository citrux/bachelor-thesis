DATE := $(shell date +%d-%m-%Y)

all: diploma.pdf

mail: diploma.pdf
	cp diploma.pdf ~/Абдрахманов_$(DATE).pdf

%.pdf: %.tex
	latexmk -pdf $<
