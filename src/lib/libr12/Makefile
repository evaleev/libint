TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

include $(TOPDIR)/src/lib/MakeVars

NAME = libr12
TARGET = $(NAME).$(LIBSUF)
COMPILER = build_$(NAME)

default: $(TOPDIR)/lib/$(TARGET)

$(TOPDIR)/lib/$(TARGET): source
	cd $(LIBSRCLINK); $(MAKE) TOPDIR=$(TOPOBJDIR) SRCDIR=$(SRCDIR)

install:: install_inc install_target

install_inc:: $(TOPDIR)/lib/$(TARGET)
	$(INSTALL) $(INSTALLDIROPT) $(includedir)/$(NAME)
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/include/$(NAME)/*.h $(includedir)/$(NAME)

install_target:: $(TOPDIR)/lib/$(TARGET)
	$(INSTALL) $(INSTALLDIROPT) $(libdir)
	$(INSTALL) $(INSTALLLIBOPT) $< $(libdir)

clean: make_libsrcdir $(LIBSRCDIR)/Makefile
	cd $(LIBSRCLINK); $(MAKE) TOPDIR=$(TOPOBJDIR) SRCDIR=$(SRCDIR) clean

oclean: make_libsrcdir $(LIBSRCDIR)/Makefile
	cd $(LIBSRCLINK); $(MAKE) TOPDIR=$(TOPOBJDIR) SRCDIR=$(SRCDIR) oclean

distclean: make_libsrcdir $(LIBSRCDIR)/Makefile
	cd $(LIBSRCLINK); $(MAKE) TOPDIR=$(TOPOBJDIR) SRCDIR=$(SRCDIR) distclean

targetclean: make_libsrcdir $(LIBSRCDIR)/Makefile
	cd $(LIBSRCLINK); $(MAKE) TOPDIR=$(TOPOBJDIR) SRCDIR=$(SRCDIR) targetclean

realclean: make_libsrcdir
	-rm -rf $(LIBSRCDIR)
	-rm -f $(LIBSRCLINK)
	-rm -f source
	-rm -f make_libsrcdir

source: make_libsrcdir $(LIBSRCDIR)/Makefile $(LIBSRCDIR)/$(NAME).h
	echo "Last time library source was generated:\c" > $@
	date >> $@

make_libsrcdir: $(TOPOBJDIR)/src/bin/$(NAME)/$(COMPILER)
	if [ -e $@ ]; then \
$(RM) -rf $(LIBSRCDIR); \
$(RM) -f $(LIBSRCLINK); \
$(RM) -f source; \
$(RM) -f make_libsrcdir; \
fi; \
if ! mkdir $(LIBSRCDIR); then \
echo Could not mkdir $(LIBSRCDIR); \
exit 1; \
fi; \
if ! $(LN_S) $(LIBSRCDIR) $(LIBSRCLINK); then \
echo Could not ln -s $(LIBSRCDIR) $(LIBSRCLINK); \
exit 1; \
fi;
	echo "Last time the library source directory was created:\c" > $@
	date >> $@

$(LIBSRCDIR)/Makefile: make_libsrcdir $(SRCDIR)/Makefile.library
	cp -f $(SRCDIR)/Makefile.library $@

$(LIBSRCDIR)/$(NAME).h: make_libsrcdir $(TOPOBJDIR)/src/bin/$(NAME)/$(COMPILER) $(LIBSRCDIR)/Makefile
	cd $(LIBSRCLINK); $(TOPOBJDIR)/src/bin/$(NAME)/$(COMPILER)
