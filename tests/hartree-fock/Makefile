TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/tests/MakeVars
-include $(TOPDIR)/src/lib/libint/MakeVars.features

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCDIR)/$(TOPDIR)/src/lib/libint -DSRCDATADIR=\"$(SRCDIR)/$(TOPDIR)/lib/basis\"

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -lint2
vpath %.a $(TOPDIR)/lib:$(TOPDIR)/lib/.libs

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M


TEST1 = hartree-fock
CXXTEST1SRC = $(TEST1).cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.$(OBJSUF))
CXXTEST1DEP = $(CXXTEST1SRC:%.cc=%.$(DEPSUF))

TEST2 = hartree-fock++
CXXTEST2SRC = $(TEST2).cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.$(OBJSUF))
CXXTEST2DEP = $(CXXTEST2SRC:%.cc=%.$(DEPSUF))

check:: check1 check2

check1::
check2::

ifeq ($(CXXGEN_SUPPORTS_CPP11),yes)
 ifeq ($(LIBINT_SUPPORTS_ONEBODY),yes)
  ifeq ($(LIBINT_SUPPORTS_ERI),yes)
   ifeq ($(LIBINT_HAS_EIGEN),yes)
    ifeq ($(LIBINT_CONTRACTED_INTS),yes)
     ifeq ($(LIBINT_SHELL_SET),1)
check1:: $(TEST1)
	./$^ $(SRCDIR)/h2o.xyz | $(PYTHON) $(SRCDIR)/$^-validate.py

check2:: $(TEST2)
	./$^ $(SRCDIR)/h2o_rotated.xyz | $(PYTHON) $(SRCDIR)/$^-validate.py ../../src/lib/libint/MakeVars.features
     endif
    endif
   endif
  endif
 endif
endif

$(TEST1): $(CXXTEST1OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS)

$(TEST2): $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) -lpthread

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST1) $(TEST2): CXX=$(CXXGEN)
$(TEST1) $(TEST2): CXXFLAGS=$(CXXGENFLAGS)
$(TEST1) $(TEST2): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST1) $(TEST2) *.o *.d

distclean:: realclean
	-rm -rf $(TOPDIR)/include/libint2/boost

realclean:: clean

targetclean:: clean

$(TOPDIR)/include/libint2/boost/preprocessor.hpp: $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

depend:: $(CXXTEST1DEP) $(CXXTEST2DEP)

ifneq ($(DODEPEND),no)
ifneq ($(CXXDEPENDSUF),none)
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

-include $(CXXTEST1DEP)
-include $(CXXTEST2DEP)
else

%.cc:: $(TOPDIR)/include/libint2/boost/preprocessor.hpp

endif
