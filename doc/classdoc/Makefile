TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

include $(TOPDIR)/doc/MakeVars
include $(TOPDIR)/doc/MakeRules
include $(TOPDIR)/doc/MakeSuffixRules

ifdef DOXYGEN
default: html/index.html

html/index.html: doxygen.cfg
	$(DOXYGEN) $^
else
default:
endif
