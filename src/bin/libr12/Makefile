TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(TOPDIR)/src/bin/MakeVars

BUILD_BINARY = build_libr12
BUILD_SRC = build_libr12.c emit_vrr_r_build.c \
emit_vrr_t1_build.c emit_vrr_t2_build.c \
emit_hrr_t_build.c mem_man.c emit_grt_order.c
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

distclean: oclean

targetclean:
	-rm -f $(BUILD_BINARY)

realclean: oclean

