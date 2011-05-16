TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(TOPDIR)/Makedirlist
include $(TOPDIR)/doc/MakeVars

ifeq ($(LIBINT_MAJOR_VERSION),2)
SUBDIRS = classdoc installation
else
SUBDIRS = progman
endif

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
	-rm -f MakeRules MakeSuffixRules MakeVars

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

