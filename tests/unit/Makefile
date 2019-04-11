TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/tests/MakeVars
-include $(TOPDIR)/src/lib/libint/MakeVars.features

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint
# include the data directory path as a definition
CPPFLAGS += -DSRCDATADIR=\"$(SRCDIR)/$(TOPDIR)/lib/basis\"
# initialize statics explicitly in test.cc
CPPFLAGS += -DLIBINT2_CONSTEXPR_STATICS=0

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -lint2
vpath %.a $(TOPDIR)/lib:$(TOPDIR)/lib/.libs

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M

TEST1 = test
CXXTEST1SRC = test.cc test-core.cc test-permute.cc \
  test-1body.cc test-boys.cc test-c-api.cc test-shell-order.cc \
  test-util.cc c-api-util.cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.$(OBJSUF))
CXXTEST1DEP = $(CXXTEST1SRC:%.cc=%.$(DEPSUF))
CTEST1SRC = c-api.c
CTEST1OBJ = $(CTEST1SRC:%.c=%.$(OBJSUF))

check:: check1

check1::

ifeq ($(CXXGEN_SUPPORTS_CPP11),yes)
 ifeq ($(LIBINT_SUPPORTS_ONEBODY),yes)
  ifeq ($(LIBINT_SUPPORTS_ERI),yes)
   ifeq ($(LIBINT_HAS_EIGEN),yes)
    ifeq ($(LIBINT_CONTRACTED_INTS),yes)
     ifeq ($(LIBINT_SHELL_SET),1)

$(TEST1): $(CXXTEST1OBJ) $(CTEST1OBJ) $(COMPUTE_LIB)
	$(LD) -o $@ $(CXXFLAGS) $(LDFLAGS) $^ $(SYSLIBS)

check1:: $(TEST1)
	./$^

     endif
    endif
   endif
  endif
 endif
endif

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST1): CXX=$(CXXGEN)
$(TEST1): CXXDEPEND=$(CXXGEN)
$(TEST1): CXXFLAGS=$(CXXGENFLAGS)
$(TEST1): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST1) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean

depend:: $(CXXTEST1DEP)

#########
# unpack bundled boost if needed
#########
ifneq ($(HAVE_SYSTEM_BOOST_PREPROCESSOR_VARIADICS),1)

$(TOPDIR)/include/libint2/boost/preprocessor.hpp: $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

$(CXXTEST1SRC):: $(TOPDIR)/include/libint2/boost/preprocessor.hpp

endif

#########
# extract header dependencies
#########
ifneq ($(DODEPEND),no)
ifneq ($(CXXDEPENDSUF),none)
%.d:: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXGENFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d:: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXGENFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

-include $(CXXTEST1DEP)
endif
