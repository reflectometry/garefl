# This program is public domain.

sinclude ../Makeconf

OBJECTS = model.o profile.o interface.o reflectivity.o resolution.o \
	  pars.o data.o fit.o fresnel.o reflrough.o

ifeq ($(ENABLE_MAGNETIC),yes)
  OBJECTS:=$(OBJECTS) r4xa.o magnetic.o

  # FIXME ENABLE_MAGNETIC is a code word for no fortran compiler available.
  # Suppress bspline and dgtsv if that is the case
  OBJECTS:=$(OBJECTS) dgtsv.o bspline.o
endif

# To crosscheck the matrix formalism against Parrat, uncomment
# the following then do:
#    cd octave; make clean; make ; octave
#    crosscheck(10)
# CXXFLAGS:=$(CXXFLAGS) -DUSE_PARRAT

TARGETS=librefl.a

.PHONY: distclean clean all

all: $(TARGETS)

librefl.a: $(OBJECTS)
	$(AR) $(ARFLAGS) librefl.a $(OBJECTS)
	$(RANLIB) librefl.a

r4xa.o: r4xa.f

dgtsv.o: dgtsv.f

$(OBJECTS): reflconfig.h

data.o: $(DATA_H)
fit.o: $(MODEL_H) $(INTERFACE_H) $(PROFILE_H) $(REFLCALC_H) $(DATA_H) \
	$(PARS_H) $(FIT_H)
pars.o: $(PARS_H)
interface.o: $(INTERFACE_H)
model.o: $(INTERFACE_H) $(MODEL_H) $(PROFILE_H)
profile.o: $(PROFILE_H)
resolution.o reflectivity.o magnetic.o: $(REFLCALC_H)
bspline.o: $(BSPLINE_H)

bspline_test$(EXE): bspline.c
	$(CC) -Wall -O3 -DTEST bspline.c -obspline_test$(EXE) -lm

check:

clean:
	-$(RM) *.o core librefl.a

distclean:
	-$(RM) *.o core *~
	-$(RM) reflconfig.h

# $Id$
