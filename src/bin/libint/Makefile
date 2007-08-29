TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/src/bin/MakeVars

BUILD_BINARY = build_libint
BUILD_SRC = build_libint.c emit_vrr_build.c emit_hrr_build.c mem_man.c \
emit_order.c emit_vrr_build_macro.c emit_hrr_build_macro.c
CONFIG = ../libint_config.h

LD = $(CC)
vpath %.a $(libdir)

BINOBJ = $(BUILD_SRC:%.c=%.o)

$(BUILD_BINARY): $(BINOBJ)
	$(LD) $^ $(LDFLAGS) $(SYSLIBS) -o $@

%.$(OBJSUF): %.c $(CONFIG)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< $(OUTPUT_OPTION)

install:

install_inc:

install_target:

uninstall:

clean: oclean
	-rm -f $(BUILD_BINARY)

oclean:
	-rm -f $(BINOBJ)

distclean: clean

targetclean:

realclean: clean
