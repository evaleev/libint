TOPDIR = ..
SRCDIR = ..

-include ../MakeSuffixRules
-include ../MakeVars
-include ../MakeVars.features

FCFLAGS := -I../include -I../include/libint2 -D__COMPILING_LIBINT2=1 $(FCFLAGS)
COMPUTE_LIB = -L../lib -lint2

.PHONY: clean distclean default make_test check_test unpack_boost

default:: fortran_example check_test libint_f.o

libint_f.o: libint2_types_f.h fortran_incldefs.h

libint2_types_f.h: $(TOPDIR)/include/libint2.h
	$(CXXCPP) $(CPPFLAGS) $< | $(PYTHON) c_to_f.py Libint_t > $@

fortran_example.o: libint_f.o

fortran_incldefs.h: $(TOPDIR)/include/libint2_types.h
	grep '^#' $< | grep -v '#include' > $@

fortran_example: fortran_example.o libint_f.o
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) $(FCLIBS) -o $@ $^

make_test::

check_test::

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
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) $(FCLIBS) -o $@ $^

test.o :: unpack_boost

unpack_boost::
	gunzip -c $(SRCDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

clean::
	-rm -f *.o *.mod fortran_example libint2_types_f.h test fortran_incldefs.h

distclean:: clean
