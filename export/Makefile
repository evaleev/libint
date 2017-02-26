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
	(cd $(TOPDIR) && $(TAR) -czf $(EXPORTDIR).tgz $(EXPORTDIR) && rm -rf $(EXPORTDIR))

$(TOPDIR)/$(EXPORTDIR)/configure:: exportdir strip_some_preproc_symb
	(cd $(TOPDIR)/$(EXPORTDIR) && aclocal -I lib/autoconf && autoconf && rm -rf autom4* aclocal.m4)

# some preprocessor symbols need to be stripped/replaced since they can be redefined by the user of the exported library
strip_some_preproc_symb:: exportdir
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h | sed -e '/^#ifndef LIBINT2_ALIGN_SIZE/,/^#endif/d' | sed -e '/^#ifndef LIBINT2_REALTYPE/,/^#endif/d' > $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.in | grep -v LIBINT_ALIGN_SIZE | grep -v HAVE_POSIX_MEMALIGN | grep -v  LIBINT_USER_DEFINED_REAL > $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.tmp \
 && printf "/* EXTRA DEFINES DETERMINED BY CONFIGURE OF THE EXPORTED LIBRARY */\n#ifndef _libint2_include_libint2_config_h_1\n#define _libint2_include_libint2_config_h_1\n#undef LIBINT2_ALIGN_SIZE\n#undef HAVE_POSIX_MEMALIGN\n#undef LIBINT2_REALTYPE\n#endif // header guard #2\n" >> $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.in

exportdir::
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/libint2.pc.in $(TOPDIR)/$(EXPORTDIR)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/INSTALL.export $(TOPDIR)/$(EXPORTDIR)/INSTALL
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/LICENSE.export $(TOPDIR)/$(EXPORTDIR)/LICENSE
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/README.md $(TOPDIR)/$(EXPORTDIR)/README.md
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/COPYING.LIB $(TOPDIR)/$(EXPORTDIR)/COPYING
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/CITATION $(TOPDIR)/$(EXPORTDIR)/CITATION
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeVars.export $(TOPDIR)/$(EXPORTDIR)/MakeVars.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/MakeRules.export $(TOPDIR)/$(EXPORTDIR)/MakeRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/MakeSuffixRules $(TOPDIR)/$(EXPORTDIR)/MakeSuffixRules.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.export $(TOPDIR)/$(EXPORTDIR)/Makefile
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/configure.export $(TOPDIR)/$(EXPORTDIR)/configure.ac
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/config.* $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/ltmain.sh $(TOPDIR)/$(EXPORTDIR)/bin
	-$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCTOPDIR)/bin/install-sh $(TOPDIR)/$(EXPORTDIR)/bin
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/libtool.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/lt*.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/autoconf/acx*.m4 $(TOPDIR)/$(EXPORTDIR)/lib/autoconf
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/basis
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/basis/* $(TOPDIR)/$(EXPORTDIR)/lib/basis	
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/chemistry
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/util
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/util/generated
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.hpp $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/basis.h.in $(TOPDIR)/$(EXPORTDIR)/include/libint2
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/chemistry/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/chemistry
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/util/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/util
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/util/generated/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/util/generated
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/include/libint2/config.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.in
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/include/libint2/cgshell_ordering.h $(TOPDIR)/$(EXPORTDIR)/include/libint2
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/util_types.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/src/lib/libint/MakeVars.features $(TOPDIR)/$(EXPORTDIR)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/lib/libint/*.h $(TOPDIR)/$(EXPORTDIR)/include
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.tests.export $(TOPDIR)/$(EXPORTDIR)/tests/Makefile
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/eri
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/eri/test.cc $(TOPDIR)/$(EXPORTDIR)/tests/eri/test.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/eri.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/eri.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/prep_libint2.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/prep_libint2.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/deriv_iter.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/deriv_iter.h
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock.cc $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock++.cc $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock++.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/h2o.xyz $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/h2o.xyz
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/h2o_rotated.xyz $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/h2o_rotated.xyz
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock-validate.py $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock-validate.py
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock++-validate.py $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock++-validate.py
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/doc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/progman.tex $(TOPDIR)/$(EXPORTDIR)/doc/progman.tex
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/*.cc $(TOPDIR)/$(EXPORTDIR)/doc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/progman/refs.bib $(TOPDIR)/$(EXPORTDIR)/doc/refs.bib
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/doc/progman/macros.tex $(TOPDIR)/$(EXPORTDIR)/doc/macros.tex
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/Libint_Logo3_alt.pdf $(TOPDIR)/$(EXPORTDIR)/doc/Libint_Logo3_alt.pdf
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/doc/Libint_Logo3_alt.eps $(TOPDIR)/$(EXPORTDIR)/doc/Libint_Logo3_alt.eps
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/Makefile.doc.export $(TOPDIR)/$(EXPORTDIR)/doc/Makefile.in
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/external
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/external/boost.tar.gz $(TOPDIR)/$(EXPORTDIR)/external/boost.tar.gz
	$(TOPOBJDIR)/config.status --config > $(TOPDIR)/$(EXPORTDIR)/compiler.config

