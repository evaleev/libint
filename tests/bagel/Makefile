TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# how to link against MPQC? provide valid SCCONFIG!
SCCONFIG = /usr/local/mpqc/current/bin/sc-config
SCCPPFLAGS := $(shell $(SCCONFIG) --cppflags)
SCLIBS := $(shell $(SCCONFIG) --libs)
SCLIBDIR  := $(shell $(SCCONFIG) --libdir)

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2
# include headers from the compiler directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint $(SCCPPFLAGS)
# include processed headers from the compiler directory
CPPFLAGS += -I$(TOPDIR)/src/bin
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint
# include bagel headers and libraries
# user must define BAGELDIR and BAGELOBJDIR (for now using files from build directory)
CPPFLAGS += -I$(BAGELDIR) -I$(BAGELBUILDDIR) -DINCLUDE_BAGEL=0
BAGELLIBS = -L$(BAGELINSTALLDIR)/lib -lbagel_parallel -lbagel_scf -lbagel_rysint -lbagel_util
BOOSTLIBS = -lboost_chrono -lboost_system

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -L$(TOPDIR)/lib -L$(TOPDIR)/lib/.libs -lint2

TEST = test
AUXSRC = $(TOPDIR)/src/bin/test_eri/eri.cc
CXXTESTSRC = $(TEST).cc $(AUXSRC)
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

$(TEST): $(CXXTESTOBJ)
	$(LD) -o $@ $(LDFLAGS) $^ $(BAGELLIBS) $(COMPILER_LIB) $(COMPUTE_LIB) -L$(SCLIBDIR) $(SCLIBS) $(SYSLIBS) $(BOOSTLIBS)

# Source files for timer and tester are to be complied using CXXGEN
$(TEST): CXX=$(CXXGEN)
$(TEST): CXXFLAGS=$(CXXGENFLAGS)
$(TEST): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST) *.o *.d

realclean:: clean

