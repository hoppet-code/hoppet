# Makefile generated automatically with
# ../scripts/makef95makefile nnpdf_hoppet -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	nnpdf_hoppet

ALLPROG = 	nnpdf_hoppet

# This will be used one day...
ALLPROGSRC =	nnpdf_hoppet.f90

ALLPROGOBJ =	nnpdf_hoppet.o

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

nnpdf_hoppet: nnpdf_hoppet.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o nnpdf_hoppet nnpdf_hoppet.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d *~ *.res fort* nnpdf_hoppet

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d 

make:
	./mkmk

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

nnpdf_hoppet.o: 
