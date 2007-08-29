TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

BUILD_BINARY = build_libint
BUILD_CSRC = build_libint.c emit_vrr_build.c emit_hrr_build.c mem_man.c \
emit_order.c emit_vrr_build_macro.c emit_hrr_build_macro.c
CONFIG = ../libint_config.h

-include $(TOPDIR)/src/bin/MakeVars
-include $(TOPDIR)/src/bin/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules

%.$(OBJSUF):: $(CONFIG)

