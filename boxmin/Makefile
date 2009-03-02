sinclude ../../reflectometry/Makeconf
CFLAGS?=-O2

%.o: %.c ; $(CC) $(CFLAGS) -c  $< -o $@

#HAVE_LAPACK = -DHAVE_LAPACK

AR?=ar
ARFLAGS?=rc
RANLIB?=ranlib
SUBDIRS=examples
.PHONY: $(SUBDIRS)

OBJS=amoeba.o nlls.o lm.o lmbc.o Axb.o misc.o
TARGET=libboxmin.a

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

amoeba.o: amoeba.c amoeba.h

nlls.o: nlls.c nlls.h
	$(CC) $(CFLAGS) -c -Ilevmar $< -o $@

lm.o: levmar/lm.c levmar/lm_core.c levmar/lm_uct.c levmar/lm.h levmar/misc.h
	$(CC) $(CFLAGS) $(HAVE_LAPACK) -c $< -o $@

lmbc.o: levmar/lmbc.c levmar/lmbc_core.c levmar/lm.h levmar/misc.h
	$(CC) $(CFLAGS) $(HAVE_LAPACK) -c $< -o $@

Axb.o: levmar/Axb.c levmar/Axb_core.c levmar/lm.h levmar/misc.h
	$(CC) $(CFLAGS) $(HAVE_LAPACK) -c $< -o $@

misc.o: levmar/misc.c levmar/lm.h levmar/misc.h
	$(CC) $(CFLAGS) $(HAVE_LAPACK) -c $< -o $@

clean: $(SUBDIRS)
	$(RM) *.o $(TARGET)

distclean: $(SUBDIRS)
	$(RM) *.o $(TARGET) *~

check:

$(SUBDIRS):
	cd $@ && $(MAKE) $(MAKECMDGOALS)

# $Id: Makefile 39 2007-08-20 19:23:20Z ziwen $
