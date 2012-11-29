TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include
# include headers from the compiler directory (should go away soon)
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint
# include generated headers from the current directory
CPPFLAGS += -I.

LIBINTDIR = $(TOPDIR)/src/bin/libint
LIBINT = $(LIBINTDIR)/libINT.a
CGSHELLOBJ = $(LIBINTDIR)/gauss.o $(LIBINTDIR)/iter.o $(LIBINTDIR)/policy.o \
$(LIBINTDIR)/policy_spec.o $(LIBINTDIR)/default_params.o

#
# making target1
#
COMP1 = generate_t0g12_code
COMP1SRC = $(COMP1).cc
COMP1OBJ = $(COMP1SRC:%.cc=%.o)

$(COMP1OBJ): CXXFLAGS := -pg $(CXXFLAGS)
$(COMP1): LDFLAGS := -pg $(LDFLAGS)

$(COMP1): $(COMP1OBJ) $(LIBINTDIR)/buildtest.o $(LIBINT)
	cd ../libint; $(MAKE) $(LIBINT)
	$(LD) -o $@ $(LDFLAGS) $^

$(COMP1OBJ):: $(LIBINTDIR)/buildtest.h

#
# making target2
#
COMP2 = profile_Gm
COMP2SRC = $(COMP2).cc
COMP2OBJ = $(COMP2SRC:%.cc=%.o)

$(COMP2): $(COMP2OBJ)
	$(CXX) -o $@ $(LDFLAGS) -lboost_chrono -lboost_system $^

$(COMP2OBJ):: $(TOPDIR)/src/lib/libint/boys.h

clean::
	-rm -rf $(COMP1) $(COMP2) *.o *.d

realclean:: clean
	-rm -rf *.h *.cc *.dot *.log

distclean:: realclean
