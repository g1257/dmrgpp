
LDFLAGS =       -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../src
CXX = g++ -pg -g3 -DNDEBUG
all: $(EXENAME)
HEADERSH = ../src/AkimaSpline.h ../src/Vector.h ../src/SampleCRSMatrix.h

all: testCRSMatrix

testCRSMatrix: $(HEADERSH) testCRSMatrix.cpp
	$(CXX) -o testCRSMatrix -I../src testCRSMatrix.cpp

../src/SampleCRSMatrix.h: ../src/SampleCRSMatrix.w
	cd ../src/; ~/software/nuweb/nuweb -v -l  -s  -d  SampleCRSMatrix.w
	
clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt


######## End of Makefile ########

