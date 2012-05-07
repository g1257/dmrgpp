CPPFLAGS = -I../src/Engine -I../src/Models/Heisenberg -I../src/Models/HeisenbergSpinOneHalf  -I../src/Models/HubbardOneBand -I../src/Models/FeAsModel \
	   -IModels/Immm -I../src/Models/FeAsBasedScExtended -IModels/ExtendedHubbard1Orb -I../../PsimagLite/src -I../../PsimagLite/src/Geometry 
CXX = g++ 

all: manual.pdf

../perlmod/DoxyDocs.pm:  ../src/dmrg.o 
	cd ../src; $(CXX) $(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep
	cd ../; doxygen

manual.tex: manual.ptex ../perlmod/DoxyDocs.pm
	perl -I../perlmod ../../PsimagLite/PsimagDoc/psimagdoc.pl < manual.ptex > manual.tex

manual.pdf: manual.tex ../README.tex
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

../README.tex:
	pod2latex -out ../README.tex  ../README.pod

clean:
	rm -f ../README.tex
	rm -f manual.pdf


