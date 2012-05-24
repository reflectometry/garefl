# This program is public domain

sinclude Makeconf
SUBDIRS=boxmin model1d src examples

REPO = git://github.com/reflectometry/garefl.git
VERSION = $(shell date +%Y.%m.%d)
DISTDIR=garefl-$(VERSION)-src
DIST=release/$(DISTDIR).tar.gz

.PHONY: $(SUBDIRS) ChangeLog

all: $(SUBDIRS)

check: $(SUBDIRS)

clean: $(SUBDIRS)

distclean: $(SUBDIRS)
	$(RM) *~ Makeconf
	-$(RM) config.cache config.log config.status lib/reflconfig.h
	-$(RM) autom4te.cache/*
	-rmdir autom4te.cache

$(SUBDIRS):
	cd $@ && $(MAKE) $(MAKECMDGOALS)

ChangeLog:
	$(SVN2CL) > /dev/null

release:
	mkdir release

configure: configure.in ax_cflags_warn_all.m4 autogen.sh
	./autogen.sh

release/garefl.html: doc/garefl.html.in
	sed -e"s,@VERSION@,$(VERSION),g" $< > $@

release/garefl-$(VERSION)-notes: RELEASE-NOTES
	cp $< $@

$(DIST): configure ChangeLog release/garefl.html release/garefl-$(VERSION)-notes
	git archive --remote $(REPO) | tar -x -C $(DISTDIR) >/dev/null
	-cp ChangeLog $(DISTDIR)
	tar cf - $(DISTDIR) | gzip > $(DIST)
	rm -rf $(DISTDIR)
REPO = git://github.com/reflectometry/garefl.git
VERSION = $(shell date +%Y.%m.%d)
DISTDIR=garefl-$(VERSION)-src
DIST=release/$(DISTDIR).tar.gz

dist: release $(DIST)


# $Id$
