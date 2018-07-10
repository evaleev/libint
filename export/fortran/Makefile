TOPDIR = ..
-include ../MakeSuffixRules
-include ../MakeVars
-include ../MakeVars.features

FCFLAGS := -I../include -I../include/libint2 -D__COMPILING_LIBINT2=1 $(FCFLAGS)
COMPUTE_LIB = -L../lib -lint2

.PHONY: clean distclean default make_test check_test

default:: fortran_example check_test libint_f.o

libint_f.o: fortran_incldefs.h

fortran_example.o: libint_f.o

fortran_incldefs.h: $(TOPDIR)/include/libint2_types.h
	grep '^#' $< | grep -v '#include' > $@; \
	awk -F'[][]' '/targets\[/ { print "#define LIBINT2_MAX_NTARGETS " $$2 }' $< >>$@

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

clean::
	-rm -f *.o *.mod fortran_example test fortran_incldefs.h

distclean:: clean
