
LDFLAGS =       -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../src
CXX = g++ -pg -g3 -DNDEBUG
all: $(EXENAME)
HEADERSH = ../src/AkimaSpline.h ../src/Vector.h ../src/SampleCRSMatrix.h ../src/Minimizer.h

all: testCRSMatrix minimizer

testCRSMatrix: $(HEADERSH) testCRSMatrix.cpp
	$(CXX) -o testCRSMatrix -I../src testCRSMatrix.cpp

minimizer: $(HEADERSH) minimizer.o
	$(CXX) -o minimizer minimizer.o -lgsl  -lgslcblas -lm

minimizer.o: $(HEADERSH) minimizer.cpp
	$(CXX) -c -I../src minimizer.cpp

../src/SampleCRSMatrix.h: ../src/SampleCRSMatrix.w
	cd ../src/; ~/software/nuweb/nuweb -v -l  -s  -d  SampleCRSMatrix.w

../src/Minimizer.h: ../src/Minimizer.w
	cd ../src/; ~/software/nuweb/nuweb -v -l  -s  -d Minimizer.w
	
clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt


######## End of Makefile ########

