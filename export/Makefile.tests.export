
TOPDIR = ..
SRCDIR = ..
-include ../MakeVars
-include ../MakeVars.features

CPPFLAGS := -I../include -I../include/libint2 -Ieri -Ihartree-fock -DPREP_LIBINT2_SKIP_BOOST -DSRCDATADIR=\"$(SRCDIR)/lib/basis\" $(CPPFLAGS)
COMPUTE_LIB = -L../lib -lint2

check:: check1 check2 check3 check4

TEST1 = eri/test
check1:: $(TEST1)
	./$(TEST1) 0 2 && ./$(TEST1) 1 1 && ./$(TEST1) 2 1

$(TEST1): eri/test.o
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) -o $@ $^

TEST2 = unit/test
TEST3 = hartree-fock/hartree-fock
TEST4 = hartree-fock/hartree-fock++

check2::
check3::
check4::

ifeq ($(LIBINT_HAS_CXX11),yes)
 ifeq ($(LIBINT_SUPPORTS_ONEBODY),yes)
  ifeq ($(LIBINT_SUPPORTS_ERI),yes)
   ifeq ($(LIBINT_HAS_EIGEN),yes)
    ifeq ($(LIBINT_CONTRACTED_INTS),yes)
     ifeq ($(LIBINT_SHELL_SET),1)
check2:: $(TEST2)
	./$^

check3:: $(TEST3)
	$^ hartree-fock/h2o.xyz | $(PYTHON) $^-validate.py

check4:: $(TEST4)
	$^ hartree-fock/h2o_rotated.xyz | $(PYTHON) $^-validate.py $(SRCDIR)/MakeVars.features
     endif
    endif
   endif
  endif
 endif
endif

CXXTEST2SRC = $(TEST2).cc $(TEST2)-core.cc $(TEST2)-permute.cc $(TEST2)-1body.cc $(TEST2)-boys.cc $(TEST2)-c-api.cc unit/c-api-util.cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.o)
CXXTEST2DEP = $(CXXTEST2SRC:%.cc=%.$(DEPSUF))
CTEST2SRC = unit/c-api.c
CTEST2OBJ = $(CTEST2SRC:%.c=%.o)


$(TEST2): $(CXXTEST2OBJ) $(CTEST2OBJ)
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) -o $@ $^

$(TEST3): $(TEST3).o
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) -o $@ $^

$(TEST4): $(TEST4).o
	$(LTLINK) $(CXX) $(LDFLAGS) $(COMPUTE_LIB) $(SYSLIBS) -lpthread -o $@ $^

.PHONY: unpack_boost

$(TEST2).o :: unpack_boost
$(TEST3).o :: unpack_boost
$(TEST4).o :: unpack_boost

unpack_boost::
	gunzip -c $(SRCDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

clean:: targetclean
	-rm -f */*.o */*.d

targetclean::
	-rm -f $(TEST1) $(TEST2) $(TEST3) $(TEST4)

distclean:: clean

