CPPFLAGS = -I../src/Engine -I../src/Models/Heisenberg -I../src/Models/HeisenbergSpinOneHalf  -I../src/Models/HubbardOneBand -I../src/Models/FeAsModel \
	   -I../src/Models/FeasBasedScExtended -I../../PsimagLite/src -I../../PsimagLite/src/Geometry
CXX = g++ 

all: manual.pdf

../perlmod/DoxyDocs.pm:  ../src/dmrg.o 
	cd ../src; $(CXX) $(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep
	cd ../; doxygen

manual.tex: manual.ptex ../perlmod/DoxyDocs.pm
	perl -I../perlmod ../../PsimagLite/PsimagDoc/psimagdoc.pl < manual.ptex > manual.tex

manual.pdf: manual.tex
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex


