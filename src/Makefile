# Makefile generated automatically with
# /ada1/lpthe/salam/scripts/makePNEW.perl libdglap.a
# default program to compile
PROG =	libdglap.a

ALLPROG = 	

# This will be used one day...
ALLPROGSRC =	

ALLPROGOBJ =	

SRCS =	hplog.f xpij2e.f xpns2e.f assertions.f90 coefficient_functions.f90 \
	convolution.f90 dglap_choices.f90 dglap_holders.f90 dglap_objects.f90 \
	evolution.f90 f77_pdftab.f90 hoppet_v1.f90 integrator.f90 \
	interpolation.f90 new_as.f90 pdf_general.f90 pdf_representation.f90 \
	pdf_tabulate.f90 qcd.f90 qcd_coupling.f90 random.f90 runge_kutta.f90 \
	sort.f90 special_functions.f90 splitting_functions.f90 \
	splitting_functions_nnlo.f90 types.f90 warnings_and_errors.f90 \
	welcome_message.f90 xa2hgp.f90 xpij2n.f90 xpij2p.f90 xpns2n.f90 \
	xpns2p.f90

POSTSRCS =	

OBJS =	hplog.o xpij2e.o xpns2e.o assertions.o coefficient_functions.o \
	convolution.o dglap_choices.o dglap_holders.o dglap_objects.o \
	evolution.o f77_pdftab.o hoppet_v1.o integrator.o interpolation.o \
	new_as.o pdf_general.o pdf_representation.o pdf_tabulate.o qcd.o \
	qcd_coupling.o random.o runge_kutta.o sort.o special_functions.o \
	splitting_functions.o splitting_functions_nnlo.o types.o \
	warnings_and_errors.o welcome_message.o xa2hgp.o xpij2n.o xpij2p.o \
	xpns2n.o xpns2p.o

POSTOBJS =	
POSTLIB = 

LIBS =  

CC = cc
CFLAGS = -O
FC = lf95
FFLAGS = -O 
F90 = lf95
F90FLAGS = -O 
LDFLAGS =  

all: libdglap.a

libdglap.a: $(OBJS)
	ar cru libdglap.a $(OBJS)
	ranlib libdglap.a

ALL: libdglap.a  $(ALLPROG)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f libdglap.a $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	/ada1/lpthe/salam/scripts/makePNEW.perl  libdglap.a

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

xpij2e.o: qcd.o
xpns2e.o: qcd.o
assertions.o: types.o
coefficient_functions.o: convolution.o qcd.o types.o
convolution.o: assertions.o integrator.o interpolation.o sort.o types.o \
	warnings_and_errors.o
dglap_holders.o: assertions.o coefficient_functions.o convolution.o \
	dglap_choices.o dglap_objects.o pdf_representation.o qcd.o types.o \
	warnings_and_errors.o
dglap_objects.o: assertions.o convolution.o dglap_choices.o \
	pdf_representation.o qcd.o splitting_functions.o types.o \
	warnings_and_errors.o
evolution.o: assertions.o convolution.o dglap_choices.o dglap_holders.o \
	dglap_objects.o pdf_representation.o qcd.o qcd_coupling.o \
	runge_kutta.o types.o warnings_and_errors.o
f77_pdftab.o: convolution.o dglap_choices.o dglap_holders.o dglap_objects.o \
	pdf_general.o pdf_tabulate.o types.o warnings_and_errors.o
hoppet_v1.o: convolution.o dglap_choices.o dglap_holders.o dglap_objects.o \
	evolution.o pdf_general.o pdf_representation.o pdf_tabulate.o \
	qcd_coupling.o types.o warnings_and_errors.o
integrator.o: types.o
interpolation.o: types.o warnings_and_errors.o
new_as.o: assertions.o qcd.o runge_kutta.o types.o warnings_and_errors.o
pdf_general.o: convolution.o pdf_representation.o types.o
pdf_representation.o: assertions.o random.o types.o warnings_and_errors.o
pdf_tabulate.o: convolution.o dglap_holders.o dglap_objects.o evolution.o \
	interpolation.o pdf_general.o pdf_representation.o qcd_coupling.o \
	types.o warnings_and_errors.o
qcd.o: types.o
qcd_coupling.o: assertions.o new_as.o types.o warnings_and_errors.o
random.o: types.o
runge_kutta.o: types.o
sort.o: assertions.o types.o warnings_and_errors.o
special_functions.o: types.o
splitting_functions.o: coefficient_functions.o convolution.o qcd.o \
	special_functions.o splitting_functions_nnlo.o types.o \
	warnings_and_errors.o
splitting_functions_nnlo.o: convolution.o dglap_choices.o qcd.o types.o \
	warnings_and_errors.o xpij2e.o xpij2n.o xpij2p.o xpns2e.o xpns2n.o \
	xpns2p.o
types.o: 
warnings_and_errors.o: types.o
