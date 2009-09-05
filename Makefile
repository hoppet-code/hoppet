
DIRS= src example_f90 benchmarking benchmarking/test_acc

all: $(DIRS)

src: force_look
	cd src; $(MAKE) $(MFLAGS)

example_f90: src force_look
	cd example_f90; $(MAKE) $(MFLAGS)

benchmarking: src force_look
	cd benchmarking; $(MAKE) $(MFLAGS)

benchmarking/test_acc: src force_look
	cd benchmarking/test_acc; $(MAKE) $(MFLAGS)

# only one directory needs installation
install: force_look
	cd src; $(MAKE) $(MFLAGS) install

check: example_f90
	scripts/check

clean: force_look
	@echo cleaning up in .
	-for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

realclean: force_look
	-for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) realclean ); done

distclean: force_look
	-for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) distclean ); done
	rm -f config.log


config.log:
	@-echo "you must run ./configure before make (try ./configure -h to see options)"
	@exit 1

force_look: config.log
	@true


