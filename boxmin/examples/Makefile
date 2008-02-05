# This program is public domain

#LAPACKLIBS = -llapack -lblas -lg2c 
CFLAGS = -Wall -O2

all: test_amoeba test_nlls

test_amoeba: test_amoeba.o ../libboxmin.a
	$(CC) $(LDFLAGS) test_amoeba.o -o $@ -L.. -lboxmin -lm

test_amoeba.o: test_amoeba.c randmtzig.ic peaks.ic ../amoeba.h
	$(CC) -c $(CFLAGS) -I.. test_amoeba.c -o test_amoeba.o

test_nlls: test_nlls.o ../libboxmin.a
	$(CC) $(LDFLAGS) test_nlls.o -o $@ -L.. -lboxmin $(LAPACKLIBS) -lm

test_nlls.o: test_nlls.c randmtzig.ic gaussfit.ic ../nlls.h
	$(CC) -c $(CFLAGS) -I.. test_nlls.c -o test_nlls.o

../libboxmin.a: builddeps ; 
builddeps: ; cd .. && $(MAKE)

clean: ; $(RM) *.o test_amoeba bowl
distclean: ; $(RM) *.o test_amoeba *~ bowl

# $Id: Makefile 7 2005-09-12 22:33:48Z pkienzle $
