TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include
# include headers from the compiler directory (should go away soon)
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include generated headers from the current directory
CPPFLAGS += -I.

ifndef COMPUTE_XX_ERI_XX_SRC
  -include Make.compute_eri
endif
CXXGENSRC = $(COMPUTE_XX_ERI_XX_SRC)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
LIBINTDIR = $(TOPDIR)/src/bin/libint
LIBINT = $(LIBINTDIR)/libINT.a
#CGSHELLOBJ = $(LIBINTDIR)/gauss.o $(LIBINTDIR)/iter.o $(LIBINTDIR)/policy.o \
#$(LIBINTDIR)/policy_spec.o $(LIBINTDIR)/default_params.o
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

ifeq ($(CXX_COMPATIBLE_WITH_CXXGEN),yes)
$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(CGSHELLOBJ)
	$(LD) -o $@ $(LDFLAGS) $^
else
$(TEST):
	@echo "Cannot test generated code! $(CXXGEN) is not compatible with $(CXXCOMP)"
endif

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ)
	$(LD) -o $@ $(LDFLAGS) $^

$(TIMEV): $(CXXTIMEVOBJ)
	$(LD) -o $@ $(LDFLAGS) $^

# Source files for timer and tester are to be complied using CXXGEN
$(TEST): CXX=$(CXXGEN)
$(TEST): CXXFLAGS=$(CXXGENFLAGS)
$(TEST): LD=$(CXXGEN)
$(TIME): CXX=$(CXXGEN)
$(TIME): CXXFLAGS=$(CXXGENFLAGS)
$(TIME): LD=$(CXXGEN)
$(TIMEV): CXX=$(CXXGEN)
$(TIMEV): CXXFLAGS=$(CXXGENFLAGS)
$(TIMEV): LD=$(CXXGEN)

$(CXXTESTOBJ): test_eri_conf.h
$(CXXTIMEOBJ): test_eri_conf.h
$(CXXERICOMPOBJ):: $(LIBINTDIR)/buildtest.h

clean::
	-rm -rf $(ERICOMP) $(TEST) $(TIME) $(TIMEV) *.o *.d

realclean:: clean
	-rm -rf Make.compute_eri *.h *.cc *.dot *.log

