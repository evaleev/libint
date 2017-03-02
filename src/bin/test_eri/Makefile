TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
-include $(TOPDIR)/src/bin/MakeVars

# Prepend this dir as include to CPPFLAGS
# (This means that the generated sources in this dir get preference over the
# global ones in $(TOPDIR)/include
$(eval CPPFLAGS = -I. $(value CPPFLAGS))
# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include
# include headers from the compiler directory (should go away soon)
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include processed headers from the compiler directory (should go away soon)
CPPFLAGS += -I$(TOPDIR)/src/bin/libint
# include headers from the library directory
CPPFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/lib/libint

ifndef COMPUTE_XX_ERI_XX_SRC
  -include Make.compute_eri
endif
CXXGENSRC = libint2_init.cc libint2_static_init.cc libint2_static_cleanup.cc $(COMPUTE_XX_ERI_XX_SRC)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
CXXGENHDRS = test_eri_conf.h libint2/libint2_types.h libint2/libint2_params.h libint2/libint2_iface.h libint2_iface_internal.h $(CXXGENSRC:%.cc=%.h)
LIBINTDIR = $(TOPDIR)/src/bin/libint
LIBINT = $(LIBINTDIR)/libINT.a
CGSHELLOBJ = $(LIBINTDIR)/gauss.o $(LIBINTDIR)/iter.o $(LIBINTDIR)/policy.o \
$(LIBINTDIR)/policy_spec.o $(LIBINTDIR)/default_params.o $(LIBINTDIR)/task.o

ERICOMP = generate_eri_code
CXXERICOMPSRC = $(ERICOMP).cc
CXXERICOMPOBJ = $(CXXERICOMPSRC:%.cc=%.o)

TEST = test_eri
CXXTESTSRC = $(TEST).cc
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

TIME = time_eri
CXXTIMESRC = $(TIME).cc
CXXTIMEOBJ = $(CXXTIMESRC:%.cc=%.o)

TEST2 = test_eri_rys
CXXTEST2SRC = $(TEST2).cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.o)

$(ERICOMP): $(CXXERICOMPOBJ) $(LIBINTDIR)/buildtest.o $(LIBINT)
	mkdir -p libint2
	cd ../libint; $(MAKE) $(LIBINT)
	$(CXX) -o $@ $^ $(SYSLIBS) $(LDFLAGS)

ifeq ($(CXX_COMPATIBLE_WITH_CXXGEN),yes)
$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(CGSHELLOBJ)
	$(CXX) -o $@ $^ $(SYSLIBS) $(LDFLAGS)
else
$(TEST):
	echo "Cannot test generated code! $(CXXGEN) is not compatible with $(CXXCOMP)"
	@exit 1
endif

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ)  $(CGSHELLOBJ)
	$(CXX) -o $@ $^ $(SYSLIBS) $(LDFLAGS)

$(TEST2): $(CXXTEST2OBJ) $(CXXGENOBJ) $(CGSHELLOBJ)
	$(CXX) -o $@ $^ $(SYSLIBS) $(LDFLAGS)

# Assume that CXXGEN has the weaker frontend of the two (there is no reason to not use it then as CXX)
# hence use it *only* to compile generated source if CXXGEN is compatible with CXX
ifeq ($(CXX_COMPATIBLE_WITH_CXXGEN),yes)
 # Only generated files to be complied using CXXGEN
 $(CXXGENOBJ): CXX=$(CXXGEN)
 $(CXXGENOBJ): CXXFLAGS=$(CXXGENFLAGS)
 $(CXXGENOBJ): LD=$(CXXGEN)
else
# Tester is not needed in this case
 #$(TEST): CXX=$(CXXGEN)
 #$(TEST): CXXFLAGS=$(CXXGENFLAGS)
 #$(TEST): LD=$(CXXGEN)
 # Source files for timer are to be complied using CXXGEN
 $(TIME): CXX=$(CXXGEN)
 $(TIME): CXXFLAGS=$(CXXGENFLAGS)
 $(TIME): LD=$(CXXGEN)
endif

$(CXXTESTOBJ): test_eri_conf.h
$(CXXTIMEOBJ): test_eri_conf.h
$(CXXERICOMPOBJ):: $(LIBINTDIR)/buildtest.h

clean::
	-rm -rf $(ERICOMP) $(TEST) $(TIME) $(TEST2) *.o *.d *.dot *.log

genclean::
	-rm -rf Make.compute_eri $(CXXGENSRC) $(CXXGENOBJ) $(CXXGENHDRS) 

realclean:: clean genclean


distclean:: realclean
	-rm -rf run_time_eri.pl run_timing_suite.pl stdtests.pl test_eri.pl time_eri.pl
