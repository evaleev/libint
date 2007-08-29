TOPDIR=../../..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

NAME = libderiv
TARGET = $(NAME).$(LIBSUF)
COMPILER = build_$(NAME)
INC = libderiv.h

-include $(TOPDIR)/src/lib/MakeVars
-include $(TOPDIR)/src/lib/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules

