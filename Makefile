DATE := $(shell date +%d-%m-%Y)
CXXFLAGS := -std=c++11

all: diploma.pdf

mail: diploma.pdf
	cp diploma.pdf ~/Абдрахманов_$(DATE).pdf

%.pdf: %.tex
	latexmk -pdf $<
