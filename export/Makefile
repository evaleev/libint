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

.PHONY: exportdir tgz strip_some_preproc_symb

export:: exportdir tgz strip_some_preproc_symb

tgz:: exportdir strip_some_preproc_symb
	(cd $(TOPDIR) && $(TAR) -czf $(EXPORTDIR).tgz $(EXPORTDIR) && rm -rf $(EXPORTDIR))

# some preprocessor symbols need to be stripped/replaced since they can be redefined by the user of the exported library
strip_some_preproc_symb:: exportdir
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h | sed -e '/^#ifndef LIBINT2_ALIGN_SIZE/,/^#endif/d' | sed -e '/^#ifndef LIBINT2_REALTYPE/,/^#endif/d' > $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h
	-rm $(TOPDIR)/$(EXPORTDIR)/include/libint2_params.h.tmp
	cat $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.in | grep -v LIBINT_ALIGN_SIZE | grep -v LIBINT_HAS_MPFR | grep -v HAVE_POSIX_MEMALIGN | grep -v  LIBINT_USER_DEFINED_REAL | grep -v "LIBINT_SHGSHELL_ORDERING " | grep -v "LIBINT_HAS_EIGEN " | grep -v "LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS " > $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.cmake.tmp \
 && printf "/* EXTRA DEFINES DETERMINED BY CONFIGURE OF THE EXPORTED LIBRARY */\n#ifndef _libint2_include_libint2_config_h_1\n#define _libint2_include_libint2_config_h_1\n#cmakedefine LIBINT2_ALIGN_SIZE @LIBINT2_ALIGN_SIZE@\n#cmakedefine HAVE_POSIX_MEMALIGN 1\n#cmakedefine LIBINT_HAS_MPFR 1\n#cmakedefine LIBINT2_REALTYPE @LIBINT2_REALTYPE@\n#cmakedefine LIBINT_SHGSHELL_ORDERING @LIBINT_SHGSHELL_ORDERING@\n#cmakedefine LIBINT_HAS_EIGEN 1\n#cmakedefine LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1\n#endif /* header guard #2 */\n" >> $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.cmake.tmp
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.cmake.tmp $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.cmake.in
	-rm $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.cmake.tmp

exportdir::
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/INSTALL.export $(TOPDIR)/$(EXPORTDIR)/INSTALL
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/LICENSE.export $(TOPDIR)/$(EXPORTDIR)/LICENSE
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/README.md $(TOPDIR)/$(EXPORTDIR)/README.md
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/INSTALL.md $(TOPDIR)/$(EXPORTDIR)/INSTALL.md
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/COPYING $(TOPDIR)/$(EXPORTDIR)/COPYING
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/COPYING.LESSER $(TOPDIR)/$(EXPORTDIR)/COPYING.LESSER
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/CITATION $(TOPDIR)/$(EXPORTDIR)/CITATION
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/basis
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/basis/* $(TOPDIR)/$(EXPORTDIR)/lib/basis
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/chemistry
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/lcao
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/util
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/include/libint2/util/generated
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2.hpp $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2
	-rm -rf $(TOPDIR)/$(EXPORTDIR)/include/libint2/basis.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/basis.h.in $(TOPDIR)/$(EXPORTDIR)/include/libint2
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/chemistry/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/chemistry
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/lcao/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/lcao
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/util/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/util
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/include/libint2/util/generated/*.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/util/generated
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/include/libint2/config.h $(TOPDIR)/$(EXPORTDIR)/include/libint2/config.h.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/libint/util_types.h $(TOPDIR)/$(EXPORTDIR)/include
	-$(INSTALL) $(INSTALLLIBOPT) $(TOPDIR)/src/lib/libint/MakeVars.features $(TOPDIR)/$(EXPORTDIR)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/lib/libint/*.h $(TOPDIR)/$(EXPORTDIR)/include
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/eri
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/eri/test.cc $(TOPDIR)/$(EXPORTDIR)/tests/eri/test.cc
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/eri.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/eri.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/prep_libint2.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/prep_libint2.h
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/unit
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/catch.hpp $(TOPDIR)/$(EXPORTDIR)/tests/unit/catch.hpp
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/fixture.h $(TOPDIR)/$(EXPORTDIR)/tests/unit/fixture.h
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/*.cc $(TOPDIR)/$(EXPORTDIR)/tests/unit
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/*.c $(TOPDIR)/$(EXPORTDIR)/tests/unit
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
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/fortran
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.F90 $(TOPDIR)/$(EXPORTDIR)/fortran
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.cc $(TOPDIR)/$(EXPORTDIR)/fortran
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/c_to_f.py $(TOPDIR)/$(EXPORTDIR)/fortran/c_to_f.py
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/make_defs.py $(TOPDIR)/$(EXPORTDIR)/fortran/make_defs.py
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/external
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/external/boost.tar.gz $(TOPDIR)/$(EXPORTDIR)/external/boost.tar.gz
	$(TOPOBJDIR)/config.status --config > $(TOPDIR)/$(EXPORTDIR)/compiler.config
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/CMakeLists.txt.export $(TOPDIR)/$(EXPORTDIR)/CMakeLists.txt
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/CMakeLists.txt.include.export $(TOPDIR)/$(EXPORTDIR)/include/CMakeLists.txt
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/cmake
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/libint2-config.cmake.in $(TOPDIR)/$(EXPORTDIR)/cmake/libint2-config.cmake.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/libint2.pc.cmake.in $(TOPDIR)/$(EXPORTDIR)/cmake/libint2.pc.cmake.in
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/hftest.cmake $(TOPDIR)/$(EXPORTDIR)/cmake/hftest.cmake
	$(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/cmake/modules
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/modules/*.cmake $(TOPDIR)/$(EXPORTDIR)/cmake/modules
	ls $(TOPDIR)/$(EXPORTDIR)/src >> $(TOPDIR)/$(EXPORTDIR)/srclist.cmake
	cp -Rap $(SRCTOPDIR)/python $(TOPDIR)/$(EXPORTDIR)/
