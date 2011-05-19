TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

TARGET = progman.$(DOCSUF)

-include $(TOPDIR)/doc/MakeVars
-include $(TOPDIR)/doc/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
