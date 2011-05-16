
TOPDIR=.
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

NAME = libint
VERSION = 2
NAMEV = $(NAME)$(VERSION)
TARGET = $(NAMEV).$(LIBSUF)

include $(TOPDIR)/MakeVars
include $(TOPDIR)/MakeRules
include $(TOPDIR)/MakeSuffixRules
VPATH = $(SRCTOPDIR)/src $(SRCTOPDIR)/include

TRUESRC = $(shell echo `find src -name '*.cc' -exec echo -n '{} ' \;`)
LIBOBJ = $(TRUESRC:%.cc=%.$(OBJSUF))

default:: $(TOPDIR)/lib/$(TARGET)

# this is how the static library is made
# NOTE: the library is made from scratch every time and the prerequisite variable is not used to avoid overflow
$(TOPDIR)/lib/$(NAMEV).a: $(LIBOBJ)
	/bin/rm -f $@
	find src -name '*.$(OBJSUF)' -print0 | xargs -0 $(AR) $(ARFLAGS) $@
	$(RANLIB) $@

# this is how shared library is made
$(TOPDIR)/lib/$(NAMEV).la: $(LIBOBJ)
	$(LTLINK) $(CXX) -o $@ $^ $(LTLINKLIBOPTS)
