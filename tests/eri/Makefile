TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2
# include headers from the compiler directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include processed headers from the compiler directory
CPPFLAGS += -I$(TOPDIR)/src/bin
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = $(TOPDIR)/lib/libint2.a

TEST = test
AUXSRC = $(TOPDIR)/src/bin/test_eri/eri.cc
CXXTESTSRC = $(TEST).cc $(AUXSRC)
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

$(TEST): $(CXXTESTOBJ)
	$(LD) -o $@ $(LDFLAGS) $^ $(COMPILER_LIB) $(COMPUTE_LIB) $(SYSLIBS)

# Source files for timer and tester are to be complied using CXXGEN
$(TEST): CXX=$(CXXGEN)
$(TEST): CXXFLAGS=$(CXXGENFLAGS)
$(TEST): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST) *.o *.d

realclean:: clean

