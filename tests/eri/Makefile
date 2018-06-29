TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/tests/MakeVars
-include $(TOPDIR)/src/lib/libint/MakeVars.features

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include
# include headers from the compiler directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include headers from the src/bin/test_eri directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/test_eri
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint
# include processed headers from the compiler directory
CPPFLAGS += -I$(TOPDIR)/src/bin/libint
# include the data directory path as a definition
CPPFLAGS += -DSRCDATADIR=\"$(SRCDIR)/$(TOPDIR)/lib/basis\"

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -lint2
vpath %.a $(TOPDIR)/lib:$(TOPDIR)/lib/.libs 

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M

TEST1 = test
CXXTEST1SRC = $(TEST1).cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.$(OBJSUF))
CXXTEST1DEP = $(CXXTEST1SRC:%.cc=%.$(DEPSUF))

TEST2 = test++
CXXTEST2SRC = $(TEST2).cc $(TEST2)-core.cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.$(OBJSUF))
CXXTEST2DEP = $(CXXTEST2SRC:%.cc=%.$(DEPSUF))

check:: check1 check2

check1:: $(TEST1)
	./$(TEST1) 0 3 && ./$(TEST1) 1 2 && ./$(TEST1) 2 1

ifeq ($(CXX_COMPATIBLE_WITH_CXXGEN),yes)
$(TEST1): $(CXXTEST1OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS)
else
$(TEST1):
	echo "Cannot test generated code! $(CXXGEN) is not compatible with $(CXXCOMP)"
	@exit 0
endif

check2::

ifeq ($(CXXGEN_SUPPORTS_CPP11),yes)
 ifeq ($(LIBINT_SUPPORTS_ONEBODY),yes)
  ifeq ($(LIBINT_SUPPORTS_ERI),yes)
   ifeq ($(LIBINT_HAS_EIGEN),yes)
    ifeq ($(LIBINT_CONTRACTED_INTS),yes)
     ifeq ($(LIBINT_SHELL_SET),1)

$(TEST2): $(CXXTEST2OBJ) $(COMPUTE_LIB)
	$(CXX) -o $@ $(LDFLAGS) $^ $(SYSLIBS)

check2:: $(TEST2)
	./$^

     endif
    endif
   endif
  endif
 endif
endif

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST1): CXX=$(CXXGEN)
$(TEST1): CXXFLAGS=$(CXXGENFLAGS)
$(TEST1): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST1) $(TEST2) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean

depend:: $(CXXTEST1DEP) $(CXXTEST2DEP)

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
