
LDFLAGS =       -lm  -lpthread /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3
CPPFLAGS = -Werror -Wall -I../src -I../src/JSON -I../src/JSON/JsonParser # -DUSE_CUSTOM_ALLOCATOR
CXX = g++ -g3
all: $(EXENAME)
HEADERSH = ../src/AkimaSpline.h ../src/Vector.h ../src/SampleCRSMatrix.h ../src/Minimizer.h
EXECUTABLES = sparseSolverTest testCRSMatrix rungeKuttaTest combineContinuedFraction \
concurrencyTest continuedFractionCollection gitrev jsonExample \
kernelPolynomial linearPrediction options randomTest range svd testLapack threads

all: $(EXECUTABLES)

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	$(CXX) $(CPPFLAGS) -c $<

testCRSMatrix: $(HEADERSH) testCRSMatrix.cpp
	$(CXX) -o testCRSMatrix -I../src testCRSMatrix.cpp

#minimizer: $(HEADERSH) minimizer.o
#	$(CXX) -o minimizer minimizer.o -lgsl  -lgslcblas -lm

#minimizer.o: $(HEADERSH) minimizer.cpp
#	$(CXX) -c -I../src minimizer.cpp

sparseSolverTest: sparseSolverTest.o
	 $(CXX) -o sparseSolverTest sparseSolverTest.o $(LDFLAGS)

../src/SampleCRSMatrix.h: ../src/SampleCRSMatrix.w
	cd ../src/; ~/software/nuweb/nuweb -v -l  -s  -d  SampleCRSMatrix.w

#../src/Minimizer.h: ../src/Minimizer.w
#	cd ../src/; ~/software/nuweb/nuweb -v -l  -s  -d Minimizer.w

clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt $(EXECUTABLES)

# dependencies brought about by Makefile.dep
 %.o: %.cpp Makefile
	$(CXX) $(CPPFLAGS) -c $<

Makefile.dep: sparseSolverTest.cpp
	$(CXX) $(CPPFLAGS) -MM sparseSolverTest.cpp  > Makefile.dep


######## End of Makefile ########

