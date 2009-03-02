include ../Makeconf

TARGETS=ga_simul.o ga.o

all: $(TARGETS)

ga_simul.o: ga_simul.c cmd_nlls.c ga.h $(MODEL1DHEADERS) $(FITHEADERS)

ga.o: ga.c ga.h

$(FITDEP):
	cd $(FIT_LIBDIR) && $(MAKE)

$(MODEL1DDEP):
	cd $(MODEL1D_LIBDIR) && $(MAKE)

VERSION ?= -$(shell date +%Y%m%d)
DISTDIR = garefl$(VERSION)
$(DISTDIR): $(GAREFLDEP) $(MODEL1DDEP) $(FITDEP)
	mkdir $(DISTDIR)
	mkdir $(DISTDIR)/bin
	cp gaplot $(DISTDIR)/bin
	mkdir $(DISTDIR)/doc
	cp doc/ga_refl.pdf $(DISTDIR)/doc
	mkdir $(DISTDIR)/lib
	cp $(GAREFLDEP) $(MODEL1DDEP) $(FITDEP) $(DISTDIR)/lib
	mkdir $(DISTDIR)/include
	cp $(MODEL1DHEADERS) $(FITHEADERS) $(DISTDIR)/include
	mkdir $(DISTDIR)/examples
	mkdir $(DISTDIR)/examples/volfrac
	mkdir $(DISTDIR)/examples/multimag
	mkdir $(DISTDIR)/examples/xray
	cp examples/volfrac/{setup.c,wc*.yor} $(DISTDIR)/examples/volfrac
	cp examples/multimag/{setup.c,f*.q[abcd]} $(DISTDIR)/examples/multimag
	cp examples/xray/{setup.c,e1085009.log} $(DISTDIR)/examples/xray
	cp examples/Makefile.dist $(DISTDIR)/examples/volfrac/Makefile
	cp examples/Makefile.dist $(DISTDIR)/examples/multimag/Makefile
	cp examples/Makefile.dist $(DISTDIR)/examples/xray/Makefile

$(DISTDIR)/Makeconf:
	echo "EXE=$(EXE)" > $(DISTDIR)/Makeconf
	echo "CC=$(CC)" >> $(DISTDIR)/Makeconf
	echo "CFLAGS=$(CFLAGS)" >> $(DISTDIR)/Makeconf
	echo "CXX=$(CXX)" >> $(DISTDIR)/Makeconf
	echo "DEFS=$(DEFS)" >> $(DISTDIR)/Makeconf
	echo "LDFLAGS=$(LDFLAGS)" >>  $(DISTDIR)/Makeconf
	echo "FLIBS=$(FLIBS)" >> $(DISTDIR)/Makeconf
	echo "LN_S=$(LN_S)" >> $(DISTDIR)/Makeconf
	echo "LIBS=$(LIBS)" >> $(DISTDIR)/Makeconf
	echo "RM=$(RM)" >> $(DISTDIR)/Makeconf

dist: $(DISTDIR) $(DISTDIR)/Makeconf
	zip -r $(DISTDIR) $(DISTDIR)
	rm -rf $(DISTDIR)

clean:
	-$(RM) *.o core $(REFLMODEL)
	rm -rf build

distclean: 
	-$(RM) *.o core *~

check:

