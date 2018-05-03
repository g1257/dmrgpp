CPPFLAGS = -I../src \
	   -I../src/Geometry
CXX = g++

all: manual.pdf

manual.tex: manual.ptex
	find ../src -iname "*.h" -or -iname "*.cpp" |\
              perl ../scripts/doc.pl manual.ptex

manual.pdf: manual.tex
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

clean:
	rm -f ../README.tex
	rm -f manual.pdf


