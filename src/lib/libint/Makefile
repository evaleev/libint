TOPDIR=../../..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

NAME = libint
TARGET = $(NAME).$(LIBSUF)
COMPILER = build_$(NAME)
INC = libint.h vrr_header.h hrr_header.h

-include $(TOPDIR)/src/lib/MakeVars
-include $(TOPDIR)/src/lib/MakeRules
-include $(TOPDIR)/src/lib/MakeSuffixRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
