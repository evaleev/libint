TOPDIR=../../..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

NAME = libint
TARGET = $(NAME).$(LIBSUF)
COMPILER = build_$(NAME)

include $(TOPDIR)/src/lib/MakeVars
include $(TOPDIR)/src/lib/MakeRules
include $(TOPDIR)/src/lib/MakeSuffixRules

