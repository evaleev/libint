TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist

SUBDIRS = src
ALLSUBDIRS = $(SUBDIRS) doc

default::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	  done

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

export::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) export) || exit 1; \
	  done
	(cd export && $(MAKE) $(DODEPENDOPT) export) || exit 1;

install::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install) || exit 1; \
	  done

install_inc::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install_inc) || exit 1; \
	  done

install_target::
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
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done

oclean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) oclean) || exit 1; \
	  done

distclean::
	for dir in $(ALLSUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) DODEPEND=no distclean) || exit 1; \
	  done
	-rm -rf autom4te.cache config.status config.log conf*.file conf*.dir *.dSYM depcheck* libtool \
Makedirlist include

targetclean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) targetclean) || exit 1; \
	  done

realclean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) realclean) || exit 1; \
	  done

check::
	for dir in tests/eri; \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) check) || exit 1; \
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
	(cd doc && $(MAKE) $(DODEPENDOPT) install) || exit 1;
