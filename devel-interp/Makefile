# Makefile generated automatically with
# ../scripts/makef95makefile test_hoppet -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	devel-interp 

ALLPROG = 	devel-interp

# This will be used one day...
ALLPROGSRC =	devel-interp.f90

ALLPROGOBJ =	devel-interp.o 	

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = g95
FFLAGS = -O -I../src 
F90 = g95
F90FLAGS = -O -I../src 
LDFLAGS =  

# Trick to enable old 'make PROG=xxx' form to still work
all: $(PROG)__

$(PROG)__: $(PROG)

ALL:  $(ALLPROG)


devel-interp: devel-interp.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o devel-interp devel-interp.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d *~ *.res fort* devel-interp  

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d 

make:
	./mkmk

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<


