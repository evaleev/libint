TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(TOPDIR)/src/bin/MakeVars

BUILD_BINARY = build_libderiv
BUILD_SRC = build_libderiv.c mem_man.c \
emit_deriv1_managers.c emit_deriv12_managers.c  emit_deriv_build.c emit_d1hrr_build.c \
emit_deriv_build_macro.c emit_d1hrr_build_macro.c
CONFIG = ../libint_config.h

LD = $(CC)
vpath %.a $(libdir)

BINOBJ = $(BUILD_SRC:%.c=%.o)

$(BUILD_BINARY): $(BINOBJ) $(SYSLIBS)
	$(LD) $^ $(LDFLAGS) -o $@

%.$(OBJSUF): %.c $(CONFIG)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< $(OUTPUT_OPTION)

install:

install_inc:

install_target:

clean: oclean

oclean:
	-rm -f $(BINOBJ)

distclean: oclean

targetclean:
	-rm -f $(BUILD_BINARY)

realclean: oclean

