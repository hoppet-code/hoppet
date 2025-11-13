


#FC = g77
FC = gfortran
FFLAGS = -O

# BELOW FOLLOW settings that need to be selected according to the f90
# compiler that was used to compile the library in ../src. If FC above
# is the same as the compiler used for the library, then F90LIB can be
# left empty (but beware, if you used a different compiler, e.g. g77,
# for compiling LHAPDF, you'll need to include a reference to the g77
# libraries (-L??? -lg2c).

# for the f90 (lf95) compiler libraries
# F90LIB= -L`which lf95 | sed s:bin/lf95:lib:` -lfj9i6 -lfj9f6 -lfj9e6 -lfccx86_6a

# for the ifort compiler libraries
#F90LIB= -L`which ifort | sed s:bin/ifort:lib:` -lifport -lifcore -limf -lm -lcxa -lirc -lunwind -lc -lirc_s

# for the g95 compiler libraries -- note that -lfrtbegin is specific
# to the combination of a g77 main program with a g95 library (ensures
# that the g77 "main"  is found before the g95 "main") -- with other compilers
# the situation may be different.
#F90LIB= -lfrtbegin -L`locate libf95.a | tail -1 | sed s:/libf95.a::` -lf95
# or the following
#F90LIB= -lfrtbegin -L$(HOME)/software/g95-install/lib/gcc-lib/i686-pc-linux-gnu/4.0.1 -lf95

# for the gfortran compiler
#F90LIB=  -L`locate libgfortran.a | tail -1 | sed s:/libgfortran.a::` -lgfortran
#F90LIB=  -L$(HOME)/software/gfortran/lib -lgfortran
#F90LIB=   -lgfortran

LIBS = -L../src -lhoppet $(F90LIB)
LHAPDFLIBS = `lhapdf-config --ldflags`

cpp_tabulation_example_qed: cpp_tabulation_example_qed.o
	$(CXX) -o cpp_tabulation_example_qed $<  $(LIBS) -lm $(F90LIB4GCC)

tabulation_example: tabulation_example.o
	$(FC) -o tabulation_example $<  $(LIBS) 

compare_lhapdf_hoppet: compare_lhapdf_hoppet.o
	$(FC) -o compare_lhapdf_hoppet $< $(LHAPDFLIBS)  $(LIBS)

convolution_example: convolution_example.o
	$(FC) -o convolution_example $< $(LHAPDFLIBS)  $(LIBS)


conveg: conveg.o
	$(FC) -o conveg $<  $(LIBS)

test_cteq_lhapdf: test_cteq_lhapdf.o
	$(FC) -o test_cteq_lhapdf $<  $(LIBS)


clean:
	rm *.o

realclean: clean
	rm convolution_example conveg  tabulation_example  cpp_tabulation_example_qed


