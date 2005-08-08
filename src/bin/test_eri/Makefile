TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CXXFLAGS += -I$(TOPDIR)/include
# include headers from the compiler directory (should go away soon)
CXXFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include generated headers from the current directory
CXXFLAGS += -I.

ifndef COMPUTE_XX_ERI_XX_SRC
  -include Make.compute_eri
endif
CXXGENSRC = $(COMPUTE_XX_ERI_XX_SRC)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
LIBINTDIR = $(TOPDIR)/src/bin/libint
LIBINT = $(LIBINTDIR)/libINT.a
CGSHELLOBJ = $(LIBINTDIR)/gauss.o $(LIBINTDIR)/iter.o $(LIBINTDIR)/policy.o \
$(LIBINTDIR)/policy_spec.o $(LIBINTDIR)/default_params.o

ERICOMP = generate_eri_code
CXXERICOMPSRC = $(ERICOMP).cc
CXXERICOMPOBJ = $(CXXERICOMPSRC:%.cc=%.o)

TEST = test_eri
AUXSRC = prep_libint2.cc eri.cc
CXXTESTSRC = $(TEST).cc $(AUXSRC)
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

TIME = time_eri
CXXTIMESRC = $(TIME).cc $(AUXSRC)
CXXTIMEOBJ = $(CXXTIMESRC:%.cc=%.o)

TIMEV = time_vecops
CXXTIMEVSRC = $(TIMEV).cc
CXXTIMEVOBJ = $(CXXTIMEVSRC:%.cc=%.o)

$(ERICOMP): $(CXXERICOMPOBJ) $(LIBINTDIR)/buildtest.o $(LIBINT)
	cd ../libint; $(MAKE) $(LIBINT)
	$(LD) -o $@ $(LDFLAGS) $^

$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(CGSHELLOBJ)
	$(CXX) -o $@ $(LDFLAGS) $^

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ) $(CGSHELLOBJ)
	$(CXX) -o $@ $(LDFLAGS) $^

$(TIMEV): $(CXXTIMEVOBJ)
	$(CXX) -o $@ $(LDFLAGS) $^

$(CXXTESTOBJ):: test_eri_conf.h
$(CXXTIMEOBJ):: test_eri_conf.h
$(CXXERICOMPOBJ):: $(LIBINTDIR)/buildtest.h

clean::
	-rm -rf $(ERICOMP) $(TEST) $(TIME) $(TIMEV) *.o *.d

realclean:: clean
	-rm -rf Make.compute_eri *.h *.cc *.dot *.log

