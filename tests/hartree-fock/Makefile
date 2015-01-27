TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCDIR)/$(TOPDIR)/src/lib/libint -DSRCDATADIR=\"$(SRCDIR)/$(TOPDIR)/lib/basis\"

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -L$(TOPDIR)/lib -L$(TOPDIR)/lib/.libs -lint2

TEST1 = hartree-fock
CXXTEST1SRC = $(TEST1).cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.o)

TEST2 = hartree-fock++
CXXTEST2SRC = $(TEST2).cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.o)

check:: check1 check2

check1:: $(TEST1)
	./$(TEST1) $(SRCDIR)/h2o.xyz | \
 awk ' {print} /^\*\*/ { eref = -74.942080057682; tol = 1e-11; if ($$5 - eref < tol && eref - $$5 < tol) {print "HF energy check passed"; exit 0} else {print "HF energy check failed"; exit 1} }'

check2::  $(TEST2)
	./$(TEST2) $(SRCDIR)/h2o.xyz | \
 awk ' {print} /^\*\*/ { eref = -76.003354058441; tol = 1e-11; if ($$5 - eref < tol && eref - $$5 < tol) {print "HF energy check passed"; exit 0} else {print "HF energy check failed"; exit 1} }'

$(TEST1): $(CXXTEST1OBJ) $(COMPILER_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(COMPUTE_LIB) $(SYSLIBS)

$(TEST2): $(CXXTEST2OBJ) $(COMPILER_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(COMPUTE_LIB) $(SYSLIBS)

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST1) $(TEST2): CXX=$(CXXGEN)
$(TEST1) $(TEST2): CXXFLAGS=$(CXXGENFLAGS)
$(TEST1) $(TEST2): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST1) $(TEST2) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean
