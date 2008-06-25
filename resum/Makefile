# Makefile generated automatically with
# ../scripts/makef95makefile tabulation_example -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	tabulation_example

ALLPROG = 	resummation_example

# This will be used one day...
ALLPROGSRC =	resummation_example.f90

ALLPROGOBJ =	resummation_example.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = ifort
FFLAGS = -O -I../src 
F90 = ifort
F90FLAGS = -O -I../src 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

resummation_example: resummation_example.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o resummation_example resummation_example.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	./mkmk

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

resummation_example.o: 
