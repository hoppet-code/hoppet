# Makefile generated automatically with
# /ada1/lpthe/salam/scripts/makePNEW.perl small_fast_tab -lcern -L../src -ldglap -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF -I../src
# default program to compile
PROG =	small_fast_tab

ALLPROG = 	lha_check_thresh small_fast_tab

# This will be used one day...
ALLPROGSRC =	lha_check_thresh.f90 small_fast_tab.f90

ALLPROGOBJ =	lha_check_thresh.o small_fast_tab.o

SRCS =	io_utils.f90 lcl_dec.f90

POSTSRCS =	

OBJS =	io_utils.o lcl_dec.o

POSTOBJS =	
POSTLIB = 

LIBS = -L/maia/cern/2000/lib -L../src -ldglap -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF  -lmathlib -lkernlib -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lg2c 

CC = cc
CFLAGS = -O
FC = ifort
FFLAGS = -g -C -I../src 
F90 = ifort
F90FLAGS = -g -C -I../src 
LDFLAGS =  

# Trick to enable old 'make PROG=xxx' form to still work
all: $(PROG)__

$(PROG)__: $(PROG)

ALL:  $(ALLPROG)

lha_check_thresh: lha_check_thresh.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o lha_check_thresh lha_check_thresh.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

small_fast_tab: small_fast_tab.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o small_fast_tab small_fast_tab.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	/ada1/lpthe/salam/scripts/makePNEW.perl  small_fast_tab -lcern -L../src -ldglap -L/ada1/lpthe/salam/utils/LHAPDF/lib -lLHAPDF -I../src

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

io_utils.o: 
lha_check_thresh.o: io_utils.o
small_fast_tab.o: io_utils.o
