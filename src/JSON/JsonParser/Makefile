#
SHELL    = /bin/bash
SRC_DIR  = .

INCLUDES =  -I. -I$(SRC_DIR) 

EXENAME		     = parseTest
CC                   = g++
cc                   = gcc
MatrixBoundsChecking = -DNDEBUG 

FLAGS        = $(MatrixBoundsChecking) -g  -Wall  -Werror  

#Note that mjson does not compile with -Werror -Wall, this is
# a temporary solution :
FLAGSC	     = $(MatrixBoundsChecking) 

OBJECTS      = main.o

all: clean showEnv $(OBJECTS) $(EXENAME)

run: all
	./parseTest small.jsn

showEnv:
	@echo ---------------------------------------------------------------- Env Info:
	@echo ------------------------------------------------ Compiler Version
	$(CC) --version	
	@echo ----------------------------------------------------------------
	@echo

main.o: main.cpp 
	$(CC) $(INCLUDES) $(FLAGS) -c main.cpp 

main.s: main.cpp 
	$(CC) -S $(INCLUDES) $(FLAGS) -c main.cpp 

$(EXENAME): $(OBJECTS)
	$(CC) $(INCLUDES) $(FLAGS) $(OBJECTS) -o $(EXENAME)  -lblas 
	cp $(EXENAME) 	$(EXENAME)_stripped
	strip $(EXENAME)_stripped

clean:
	rm -f $(EXENAME) *.o

########################End of Makefile #################################


