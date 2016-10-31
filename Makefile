
all: manuscript.pdf

manuscript.tex: manuscript.Rnw
	/bin/echo -e 'options(encoding="UTF-8")\nlibrary(knitr)\nknit("manuscript.Rnw", encoding="utf8")' | R --no-restore --no-save

manuscript.pdf: manuscript.tex manuscript.bib
	pdflatex manuscript
	biber manuscript
	pdflatex manuscript
	pdflatex manuscript
