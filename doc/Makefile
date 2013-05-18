TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist
-include $(TOPDIR)/doc/MakeVars

SUBDIRS = classdoc progman

default::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(JOBS)) || exit 1; \
	  done

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

install::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install) || exit 1; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install_target) || exit 1; \
	  done

install-pdf:: pdf
	(cd progman && make install-pdf) || exit 1

install-ps:: ps
	(cd progman && make install-ps) || exit 1

install-dvi:: dvi
	(cd progman && make install-dvi) || exit 1

install-html:: html
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) install-html) || exit 1; \
	  done

pdf::
	(cd progman && make pdf) || exit 1

ps::
	(cd progman && make ps) || exit 1

dvi::
	(cd progman && make dvi) || exit 1

html::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) html) || exit 1; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done

distclean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) distclean) || exit 1; \
	  done
	-rm -f MakeRules MakeVars

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

