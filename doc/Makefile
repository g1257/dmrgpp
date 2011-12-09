CPPFLAGS = -IEngine -IModels/Heisenberg -IModels/HeisenbergSpinOneHalf  -I../../PsimagLite/src -I../../PsimagLite/src/Geometry
CXX = g++ 

all: manual.pdf

../perlmod/DoxyDocs.pm: ../src/dmrg.o
	cd ../src; $(CXX) $(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep
	cd ../; doxygen

manual.tex: manual.ptex ../perlmod/DoxyDocs.pm
	perl -I../perlmod manual.pl < manual.ptex > manual.tex

manual.pdf: manual.tex
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex


