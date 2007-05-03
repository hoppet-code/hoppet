# Makefile generated automatically with
# ../scripts/makef95makefile tabulation_example -L../src -lhoppet_v1 -I../src -remake ./mkmk
# default program to compile
PROG =	tabulation_example

ALLPROG = 	lha_check_thresh tabulation_example

# This will be used one day...
ALLPROGSRC =	lha_check_thresh.f90 tabulation_example.f90

ALLPROGOBJ =	lha_check_thresh.o tabulation_example.o

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

lha_check_thresh: lha_check_thresh.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o lha_check_thresh lha_check_thresh.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example: tabulation_example.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example tabulation_example.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

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

lha_check_thresh.o: 
tabulation_example.o: 
