
CXX = g++



# We're using ansi C++
CPPFLAGS += -pedantic -std=c++11

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Wendif-labels

STRIP_COMMAND = gdb-add-index
CPPFLAGS +=   -g3

CPPFLAGS += -Werror



all: MuFunction

MuFunction.o: MuFunction.cpp  Makefile.mu
	$(CXX) $(CPPFLAGS) -c MuFunction.cpp

MuFunction: MuFunction.o  Makefile.mu
	$(CXX) $(LDFLAGS) -o MuFunction MuFunction.o

clean::
	rm -f MuFunction.o  MuFunction

