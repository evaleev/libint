TOPDIR=..
TO_TOPDIR=$(TOPDIR)
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))
TAR = tar

-include $(TOPDIR)/src/lib/MakeVars

default::

install::

clean::

realclean::

.PHONY: exportdir tgz

export:: exportdir $(TOPDIR)/$(EXPORTDIR)/configure tgz

tgz:: exportdir $(TOPDIR)/$(EXPORTDIR)/configure
	(cd $(TOPDIR) && $(TAR) -cvzf $(EXPORTDIR).tgz $(EXPORTDIR) && rm -rf $(EXPORTDIR))

$(TOPDIR)/$(EXPORTDIR)/configure:: exportdir strip_some_preproc_symb
	(cd $(TOPDIR)/$(EXPORTDIR) && aclocal -I lib/autoconf && autoconf && rm -rf autom4* aclocal.m4)

# some preprocessor symbols need to be stripped/replaced since they can be redefined by the user of the exported library
strip_some_preproc_symb:: exportdir
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h | sed -e '/^#ifndef LIBINT2_ALIGN_SIZE/,/^#endif/d' | sed -e '/^#ifndef LIBINT2_REALTYPE/,/^#endif/d' > $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.in | grep -v LIBINT_ALIGN_SIZE | grep -v HAVE_POSIX_MEMALIGN | grep -v  LIBINT_USER_DEFINED_REAL > $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.tmp \
 && printf "/* EXTRA DEFINES DETERMINED BY CONFIGURE OF THE EXPORTED LIBRARY */\n#ifndef _libint2_include_libint2config_h_1\n#define _libint2_include_libint2config_h_1\n#undef LIBINT2_ALIGN_SIZE\n#undef HAVE_POSIX_MEMALIGN\n#undef LIBINT2_REALTYPE\n#endif // header guard #2\n" >> $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.in

exportdir::
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/config.* $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/ltmain.sh $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/install-sh $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeVars.export $(TOPDIR)/$(EXPORTDIR)/MakeVars.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeRules.export $(TOPDIR)/$(EXPORTDIR)/MakeRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeSuffixRules.export $(TOPDIR)/$(EXPORTDIR)/MakeSuffixRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.export $(TOPDIR)/$(EXPORTDIR)/Makefile
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/configure.export $(TOPDIR)/$(EXPORTDIR)/configure.ac
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/libtool.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/lt*.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2_intrinsic_types.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2_intrinsic_operations.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/util_types.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/cgshellinfo.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/lib/libint/*.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/libint2.pc.in $(TOPDIR)/$(EXPORTDIR)
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/src/bin/libint/cgshell_ordering.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/include/libint2_config.h $(TOPDIR)/$(EXPORTDIR)/include/libint2_config.h.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/INSTALL.export $(TOPDIR)/$(EXPORTDIR)/INSTALL
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/LICENSE $(TOPDIR)/$(EXPORTDIR)/LICENSE
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/README $(TOPDIR)/$(EXPORTDIR)/README
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/COPYING $(TOPDIR)/$(EXPORTDIR)/COPYING
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/CITATION $(TOPDIR)/$(EXPORTDIR)/CITATION
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/eri/test.cc $(TOPDIR)/$(EXPORTDIR)/tests/test_eri.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/eri.h $(TOPDIR)/$(EXPORTDIR)/tests/eri.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/eri.cc $(TOPDIR)/$(EXPORTDIR)/tests/eri.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/prep_libint2.h $(TOPDIR)/$(EXPORTDIR)/tests/prep_libint2.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/deriv_iter.h $(TOPDIR)/$(EXPORTDIR)/tests/deriv_iter.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.tests.export $(TOPDIR)/$(EXPORTDIR)/tests/Makefile
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/doc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/progman.tex $(TOPDIR)/$(EXPORTDIR)/doc/progman.tex
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/*.cc $(TOPDIR)/$(EXPORTDIR)/doc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/refs.bib $(TOPDIR)/$(EXPORTDIR)/doc/refs.bib
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/doc/progman/macros.tex $(TOPDIR)/$(EXPORTDIR)/doc/macros.tex
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/Libint_Logo3_alt.pdf $(TOPDIR)/$(EXPORTDIR)/doc/Libint_Logo3_alt.pdf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/Libint_Logo3_alt.eps $(TOPDIR)/$(EXPORTDIR)/doc/Libint_Logo3_alt.eps
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.doc.export $(TOPDIR)/$(EXPORTDIR)/doc/Makefile.in
	$(TOPOBJDIR)/config.status --config > $(TOPDIR)/$(EXPORTDIR)/compiler.config

