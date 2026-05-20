FC = gfortran
FFLAGS ?= -O3 -ffixed-line-length-none
AR ?= ar
ARFLAGS ?= rcs

SOURCES := $(wildcard *.f)
OBJECTS := $(SOURCES:.f=.o)
LIB := libhplog6.a

.PHONY: all clean

all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(LIB)
