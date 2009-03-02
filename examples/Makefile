
SUBMAKEDIRS = $(dir $(wildcard */Makefile))

#$(warning SUBMAKEDIRS=$(SUBMAKEDIRS))

.PHONY: $(SUBMAKEDIRS)

all: $(SUBMAKEDIRS)

clean: $(SUBMAKEDIRS)

check: volfrac/

distclean: $(SUBMAKEDIRS)
	-$(RM) *~ 

$(SUBMAKEDIRS): ; cd $@ && $(MAKE) -k $(MAKECMDGOALS)
