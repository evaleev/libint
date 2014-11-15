TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2
# include headers from the compiler directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include headers from the src/bin/test_eri directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/test_eri
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint
# include processed headers from the compiler directory
CPPFLAGS += -I$(TOPDIR)/src/bin/libint

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -L$(TOPDIR)/lib -L$(TOPDIR)/lib/.libs -lint2

TEST = test
AUXSRC = $(TOPDIR)/src/bin/test_eri/eri.cc
CXXTESTSRC = $(TEST).cc $(AUXSRC)
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

check:: $(TEST)
	./$(TEST) 0 3 && ./$(TEST) 1 2 && ./$(TEST) 2 1

ifeq ($(CXX_COMPATIBLE_WITH_CXXGEN),yes)
$(TEST): $(CXXTESTOBJ) $(COMPILER_LIB)
	$(LD) -o $@ $(LDFLAGS) $^  $(COMPUTE_LIB) $(SYSLIBS)
else
$(TEST):
	echo "Cannot test generated code! $(CXXGEN) is not compatible with $(CXXCOMP)"
	@exit 0
endif

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST): CXX=$(CXXGEN)
$(TEST): CXXFLAGS=$(CXXGENFLAGS)
$(TEST): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean
