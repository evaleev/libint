TOPDIR = ..
SRCDIR = ..

-include ../MakeSuffixRules
-include ../MakeVars
-include ../MakeVars.features

FCFLAGS := -I../include -I../include/libint2 -D__COMPILING_LIBINT2=1 $(FCFLAGS)
COMPUTE_LIB = -L../lib -lint2
CPPFLAGS += -DPREP_LIBINT2_SKIP_BOOST

.PHONY: clean distclean default make_test check_test

default::
make_test::
check_test::

ifeq ($(ENABLE_FORTRAN),yes)

default:: fortran_example check_test

libint_f.o: libint2_types_f.h fortran_incldefs.h

libint2_types_f.h: $(TOPDIR)/include/libint2.h
	$(CXXCPP) $(CPPFLAGS) $(CXXFLAGS) $< | $(PYTHON) c_to_f.py Libint_t > $@

fortran_example.o: libint_f.o

fortran_incldefs.h: $(TOPDIR)/include/libint2_types.h
	grep '^#' $< | grep -v '#include' > $@

fortran_example: fortran_example.o libint_f.o
	$(LTLINK) $(CXX) $(CXXFLAGS) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) $(FCLIBS) -o $@ $^

ifeq ($(LIBINT_HAS_CXX11),yes)
 ifeq ($(LIBINT_SUPPORTS_ONEBODY),yes)
  ifeq ($(LIBINT_SUPPORTS_ERI),yes)
   ifeq ($(LIBINT_HAS_EIGEN),yes)
    ifeq ($(LIBINT_CONTRACTED_INTS),yes)
     ifeq ($(LIBINT_SHELL_SET),1)

make_test:: test

check_test:: make_test
	./test

     endif
    endif
   endif
  endif
 endif
endif

test: test.o test-eri.o libint_f.o
	$(LTLINK) $(CXX) $(CXXFLAGS) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) $(FCLIBS) -o $@ $^

ifneq ($(HAVE_SYSTEM_BOOST_PREPROCESSOR_VARIADICS),1)

$(TOPDIR)/include/libint2/boost/preprocessor.hpp:: $(SRCDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

test.cc test-eri.cc :: $(TOPDIR)/include/libint2/boost/preprocessor.hpp

endif

endif # ENABLE_FORTRAN=yes

clean::
	-rm -f *.o *.mod fortran_example libint2_types_f.h test fortran_incldefs.h

distclean:: clean
