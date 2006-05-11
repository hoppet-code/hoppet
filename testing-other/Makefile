# Makefile generated automatically with
# /ada1/lpthe/salam/scripts/makePNEW.perl testcteq6 -lcern -L../src -lpdfconv -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF -I../src
# default program to compile
PROG =	testcteq6

ALLPROG = 	testcteq6-v2 testcteq6

# This will be used one day...
ALLPROGSRC =	testcteq6-v2.f90 testcteq6.f90

ALLPROGOBJ =	testcteq6-v2.o testcteq6.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L/maia/cern/2000/lib -L../src -lpdfconv -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF  -lmathlib -lkernlib -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lg2c -lm 

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -O3 -march=i686 -I../src 
F90 = gfortran
F90FLAGS = -O3 -march=i686 -I../src 
LDFLAGS =  

# Trick to enable old 'make PROG=xxx' form to still work
all: $(PROG)__

$(PROG)__: $(PROG)

ALL:  $(ALLPROG)

testcteq6-v2: testcteq6-v2.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o testcteq6-v2 testcteq6-v2.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

testcteq6: testcteq6.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o testcteq6 testcteq6.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	/ada1/lpthe/salam/scripts/makePNEW.perl  testcteq6 -lcern -L../src -lpdfconv -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF -I../src

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

testcteq6-v2.o: 
testcteq6.o: 
