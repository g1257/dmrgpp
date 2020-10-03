CPPFLAGS = -I../src \
	   -I../src/Geometry
CXX = g++

all: manual.pdf tutorial.pdf

manual.tex: manual.ptex
	find ../src -iname "*.h" -or -iname "*.cpp" |\
              perl ../scripts/doc.pl manual.ptex

tutorial.tex: tutorial.ptex ../src/*h ../drivers/*cpp
	find ../src ../drivers -iname "*.h" -or -iname "*.cpp" |\
perl ../scripts/doc.pl tutorial.ptex

manual.pdf: manual.tex
	pdflatex manual.tex
#	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

tutorial.pdf: tutorial.tex
	pdflatex tutorial.tex
#	bibtex tutorial.tex
	pdflatex tutorial.tex
	pdflatex tutorial.tex

clean:
	rm -f ../README.tex
	rm -f manual.pdf manual.tex tutorial.pdf tutorial.tex


