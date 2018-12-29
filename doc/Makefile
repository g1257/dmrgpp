CPPFLAGS = -I../src/Engine -I../src/Models/Heisenberg -I../src/Models/Heisenberg -I../src/Models/HubbardOneBand -I../src/Models/FeAsModel \
	   -IModels/Immm -I../src/Models/FeAsBasedScExtended -IModels/ExtendedHubbard1Orb -I../../PsimagLite/src -I../../PsimagLite/src/Geometry 
CXX = g++ 

all: manual.pdf

FILES = $(shell find ../src -iname "*.h" -or -iname "*.cpp")
manual.tex: manual.ptex ../README.md  $(FILES)
	find ../src -iname "*.h" -or -iname "*.cpp" |\
              perl ../../PsimagLite/scripts/doc.pl manual.ptex

manual.pdf: manual.tex ../README.md ../src/*h ../src/Engine/*h ../src/*cpp
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

clean:
	rm -f ../README.tex
	rm -f manual.tex manual.pdf


