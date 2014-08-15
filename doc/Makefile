all:
	pdflatex paper.tex
	bibtex paper.aux
	pdflatex paper.tex
	pdflatex paper.tex

clean:
	$(RM) *.aux *.log *.pdf *.blg *.bbl

