# Makefile generated automatically with
# ../scripts/makef95makefile small_fast_tab -L../src -lhoppet_v1 -I../src -g77libs -L/home/salam/utils/LHAPDF/lib -lLHAPDF -remake ./mkmk
# default program to compile
PROG =	small_fast_tab

ALLPROG = 	benchmark small_fast_tab tabulation_example

# This will be used one day...
ALLPROGSRC =	benchmark.f90 small_fast_tab.f90 tabulation_example.f90

ALLPROGOBJ =	benchmark.o small_fast_tab.o tabulation_example.o

SRCS =	io_utils.f90 lcl_dec.f90 NameSelect.f90

POSTSRCS =	

OBJS =	io_utils.o lcl_dec.o NameSelect.o

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1 -L/home/salam/utils/LHAPDF/lib -lLHAPDF  -lfj9i6  -lg2c -lm 

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

benchmark: benchmark.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o benchmark benchmark.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

small_fast_tab: small_fast_tab.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o small_fast_tab small_fast_tab.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

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

benchmark.o: NameSelect.o io_utils.o
io_utils.o: 
NameSelect.o: io_utils.o
small_fast_tab.o: io_utils.o
tabulation_example.o: 
