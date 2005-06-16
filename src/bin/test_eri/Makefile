TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(TOPDIR)/src/bin/MakeVars

# include headers from the compiler directory (should go away soon)
CXXFLAGS += -I$(SRCDIR)/$(TOPDIR)/src/bin/libint
# include generated headers from the current directory
CXXFLAGS += -I.

CXXGENSRC = $(shell echo `ls _aB*.cc`)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
LIBINT = $(TOPDIR)/src/bin/libint/libINT.a

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

$(ERICOMP): $(CXXERICOMPOBJ) $(TOPDIR)/src/bin/libint/buildtest.o $(LIBINT)
	$(LD) -o $@ $^

$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

clean::
	-rm -rf $(ERICOMP) $(TEST) $(TIME) $(ERICOMPOBJ) $(CXXTESTOBJ) $(CXXTIMEOBJ) $(CXXGENOBJ)

