# Makefile generated automatically with
# ../scripts/makef95makefile tabulation_example -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	tabulation_example

ALLPROG = 	tabulation_example

# This will be used one day...
ALLPROGSRC =	tabulation_example.f90

ALLPROGOBJ =	tabulation_example.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = lf95
FFLAGS = -O -I../src 
F90 = lf95
F90FLAGS = -O -I../src 
LDFLAGS =  

# Trick to enable old 'make PROG=xxx' form to still work
all: $(PROG)__

$(PROG)__: $(PROG)

ALL:  $(ALLPROG)

tabulation_example: tabulation_example.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example tabulation_example.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	./mkmk

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

tabulation_example.o: 
