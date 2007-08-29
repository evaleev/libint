TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist

SUBDIRS = libint
ifeq ($(HAVE_LIBINT_SRC_BIN_LIBR12),yes)
SUBDIRS := $(SUBDIRS) libr12
endif
ifeq ($(HAVE_LIBINT_SRC_BIN_LIBDERIV),yes)
SUBDIRS := $(SUBDIRS) libderiv
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
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) distclean) || exit 1; \
	  done
	-rm -rf MakeVars libint_config.h

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

