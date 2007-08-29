TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

BUILD_BINARY = build_libr12
BUILD_CSRC = build_libr12.c emit_vrr_r_build.c \
emit_vrr_t1_build.c emit_vrr_t2_build.c \
emit_hrr_t_build.c mem_man.c emit_grt_order.c
CONFIG = ../libint_config.h

-include $(TOPDIR)/src/bin/MakeVars
-include $(TOPDIR)/src/bin/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules

%.$(OBJSUF):: $(CONFIG)
