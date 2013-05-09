FIGURES=fret-setup.pdf

all : project.pdf

%.tex : %.mkd
	pandoc -t latex -s $< -o $@ --bibliography=${HOME}/lori/papers/papers.bib/library.bib

%.pdf : %.tex ${FIGURES}
	pdflatex $<

%.pdf : %.svg
	inkscape --export-pdf=$@ $<