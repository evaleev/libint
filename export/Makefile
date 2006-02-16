TOPDIR=..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

include $(TOPDIR)/src/lib/MakeVars

default::

install::

clean::

realclean::

.PHONY: exportdir tgz

export:: exportdir $(TOPDIR)/$(EXPORTDIR)/configure tgz

tgz:: exportdir $(TOPDIR)/$(EXPORTDIR)/configure
	(cd $(TOPDIR) && tar -cvzf $(EXPORTDIR).tgz $(EXPORTDIR) && rm -rf $(EXPORTDIR))

$(TOPDIR)/$(EXPORTDIR)/configure:: exportdir
	(cd $(TOPDIR)/$(EXPORTDIR) && aclocal -I lib/autoconf && autoconf)

exportdir::
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/config.* $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/ltmain.sh $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/install-sh $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeVars.export $(TOPDIR)/$(EXPORTDIR)/MakeVars.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeRules.export $(TOPDIR)/$(EXPORTDIR)/MakeRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeSuffixRules.export $(TOPDIR)/$(EXPORTDIR)/MakeSuffixRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.export $(TOPDIR)/$(EXPORTDIR)/Makefile
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/configure.export $(TOPDIR)/$(EXPORTDIR)/configure.in
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/libtool.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2_types.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/libint2_config.export $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/INSTALL.export $(TOPDIR)/$(EXPORTDIR)/INSTALL
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/LICENSE $(TOPDIR)/$(EXPORTDIR)/LICENSE

