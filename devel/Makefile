# Makefile generated automatically with
# ../scripts/makef95makefile test_hoppet -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	test_hoppet test2_hoppet test_wuki test_evolution

ALLPROG = 	test_hoppet test2_hoppet test_wuki test_evolution

# This will be used one day...
ALLPROGSRC =	test_hoppet.f90 test2_hoppet.f90 test_wuki.f90 test_evolution.f90

ALLPROGOBJ =	test_hoppet.o 	test2_hoppet.o test_wuki.o test_evolution.o

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

# Trick to enable old 'make PROG=xxx' form to still work
all: $(PROG)__

$(PROG)__: $(PROG)

ALL:  $(ALLPROG)

test_evolution: test_evolution.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o test_evolution test_evolution.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

test_wuki: test_wuki.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o test_wuki test_wuki.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

test2_hoppet: test2_hoppet.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o test2_hoppet test2_hoppet.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

test_hoppet: test_hoppet.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o test_hoppet test_hoppet.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d *~ *.res fort* test_hoppet  test2_hoppet test_wuki test_evolution

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d 

make:
	./mkmk

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

test_hoppet.o: 
