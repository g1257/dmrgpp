
LDFLAGS =       -lm  -lpthread -lblas -llapack
CPPFLAGS = -Werror -Wall -I../src -I../src/JSON -I../src/JSON/JsonParser
CXX = g++ -g3 -DNDEBUG
all: $(EXENAME)
HEADERSH = ../src/AkimaSpline.h ../src/Vector.h ../src/SampleCRSMatrix.h ../src/Minimizer.h

all: testCRSMatrix minimizer

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	$(CXX) $(CPPFLAGS) -c $< 

Makefile.dep: dmrg.cpp
	$(CXX) $(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep

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

