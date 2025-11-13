
# -lgfortran not needed for the gfortran compilations, but it is needed
# for C++ compilations.
F90LIB4GCC= -L/opt/homebrew/opt/gcc/lib/gcc/current -lgfortran
#F90LIB4GCC=  -lgfortran

CXXFLAGS= -O3 -g -Wall -std=c++14 -I../src
LIBS=-L../src -lhoppet

all: tabulation_example tabulation_example_qed structure_functions_example

tabulation_example: tabulation_example.o
	$(CXX) -o tabulation_example $<  $(LIBS) -lm $(F90LIB4GCC)

tabulation_example_qed: tabulation_example_qed.o
	$(CXX) -o tabulation_example_qed $<  $(LIBS) -lm $(F90LIB4GCC)

structure_functions_example: structure_functions_example.o
	$(CXX) -o structure_functions_example $<  $(LIBS) -lm $(F90LIB4GCC)

clean:
	rm *.o
