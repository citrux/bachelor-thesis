all: bachelor-thesis.pdf

%.pdf: %.tex plots
	latexmk -pdf $<

plots: solver.py
	python solver.py
