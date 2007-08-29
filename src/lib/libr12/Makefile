TOPDIR=../../..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

NAME = libr12
TARGET = $(NAME).$(LIBSUF)
COMPILER = build_$(NAME)
INC = libr12.h

include $(TOPDIR)/src/lib/MakeVars
include $(TOPDIR)/src/lib/MakeRules
include $(TOPDIR)/src/lib/MakeSuffixRules

