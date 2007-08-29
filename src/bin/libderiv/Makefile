TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

BUILD_BINARY = build_libderiv
BUILD_CSRC = build_libderiv.c mem_man.c \
emit_deriv1_managers.c emit_deriv12_managers.c  emit_deriv_build.c emit_d1hrr_build.c \
emit_deriv_build_macro.c emit_d1hrr_build_macro.c
CONFIG = ../libint_config.h

-include $(TOPDIR)/src/bin/MakeVars
-include $(TOPDIR)/src/bin/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules

%.$(OBJSUF):: $(CONFIG)
