EXE =	virtualTumour

SRC = \
	initTumour.cpp \
	evolveTumour.cpp \
	particleCell_main.cpp

INC = \
	initTumour.hpp \
	evolveTumour.hpp

OBJ = \
	initTumour.o \
	evolveTumour.o \
	particleCell_main.o	

# if compiling for serial or openmp code
# ... on Macbook
#CXX = /usr/local/bin/g++-8
# ... on CAMP
CXX = g++  
# if compiling MPI c++ code for single machine (mpiexec from mpich) 
#CXX = mpic++

DEBUG = -g
CFLAGS = -fopenmp -O3 -std=c++11 $(DEBUG)
LFLAGS = -fopenmp

all : $(EXE)

$(EXE) : $(OBJ)
	$(CXX) -o $@  $(LFLAGS) $(OBJ)

#$(OBJ) : $(INC)
#	$(CXX) $(CFLAGS) -c $(SRC)

initTumour.o : initTumour.hpp initTumour.cpp
	$(CXX) $(CFLAGS) -c initTumour.cpp	
evolveTumour.o : initTumour.hpp evolveTumour.hpp evolveTumour.cpp
	$(CXX) $(CFLAGS) -c evolveTumour.cpp	
particleCell_main.o : initTumour.hpp evolveTumour.hpp particleCell_main.cpp
	$(CXX) $(CFLAGS) -c particleCell_main.cpp

clean:
	\rm $(OBJ) $(EXE)
