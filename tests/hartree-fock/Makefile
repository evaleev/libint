TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCDIR)/$(TOPDIR)/src/lib/libint -DSRCDATADIR=\"$(SRCDIR)/$(TOPDIR)/lib/data\"

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -L$(TOPDIR)/lib -L$(TOPDIR)/lib/.libs -lint2

TEST = hartree-fock
CXXTESTSRC = $(TEST).cc
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

check:: $(TEST)
	./$(TEST) $(SRCDIR)/h2o.xyz | \
 awk ' {print} /^\*\*/ { eref = -74.942080057682; tol = 1e-11; if ($$5 - eref < tol && eref - $$5 < tol) {print "HF energy check passed"; exit 0} else {print "HF energy check failed"; exit 1} }'

$(TEST): $(CXXTESTOBJ) $(COMPILER_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(COMPUTE_LIB) $(SYSLIBS)

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST): CXX=$(CXXGEN)
$(TEST): CXXFLAGS=$(CXXGENFLAGS)
$(TEST): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean
