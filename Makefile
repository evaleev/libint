TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist
-include $(TOPDIR)/src/lib/MakeVars

SUBDIRS = src
CHECKSUBDIRS = tests/eri tests/hartree-fock
CLEANSUBDIRS = $(SUBDIRS) $(CHECKSUBDIRS)
ALLSUBDIRS = $(CLEANSUBDIRS) doc $(CHECKSUBDIRS)

default::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	  done

all:: default

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

export::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) export) || exit 1; \
	  done
	(cd export && $(MAKE) $(DODEPENDOPT) export) || exit 1;

install:: all install_pkgconfig install_inc install_data
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install) || exit 1; \
	  done

ifdef pkgconfigdir
install_pkgconfig:: 
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(pkgconfigdir)
	$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/libint2.pc $(DESTDIR)$(pkgconfigdir)
endif

install_data::
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(datadir)/basis
	$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/basis/* $(DESTDIR)$(datadir)/basis

install_inc:: all
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install_inc) || exit 1; \
	  done

install_target:: all
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install_target) || exit 1; \
	  done

uninstall::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) uninstall) || exit 1; \
	  done

clean::
	for dir in $(CLEANSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done

oclean::
	for dir in $(CLEANSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) oclean) || exit 1; \
	  done

distclean::
	for dir in $(ALLSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) distclean) || exit 1; \
	  done
	-rm -rf autom4te.cache config.status config.log conf*.file conf*.dir *.dSYM depcheck* libtool \
Makedirlist libint2.pc include/libint2
	test -f libint2.pc.in || rm -rf include/libint2

targetclean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) targetclean) || exit 1; \
	  done

realclean::
	for dir in $(CLEANSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) realclean) || exit 1; \
	  done

check::
	for dir in $(CHECKSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) check) || exit 1; \
	  done

install-pdf:: pdf
	(cd doc && $(MAKE) $(DODEPENDOPT) install-pdf) || exit 1;

install-ps:: ps
	(cd doc && $(MAKE) $(DODEPENDOPT) install-ps) || exit 1;

install-dvi:: dvi
	(cd doc && $(MAKE) $(DODEPENDOPT) install-dvi) || exit 1;

install-html:: html
	(cd doc && $(MAKE) $(DODEPENDOPT) install-html) || exit 1;

pdf::
	(cd doc && $(MAKE) $(DODEPENDOPT) pdf) || exit 1;

ps::
	(cd doc && $(MAKE) $(DODEPENDOPT) ps) || exit 1;

dvi::
	(cd doc && $(MAKE) $(DODEPENDOPT) dvi) || exit 1;

html::
	(cd doc && $(MAKE) $(DODEPENDOPT) html) || exit 1;

doc::
	(cd doc && $(MAKE) $(DODEPENDOPT) install-pdf install-html) || exit 1;
