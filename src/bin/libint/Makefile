
# Define LOCAL to compile directly in the source directory
ifndef LOCAL
  TOPDIR=../../..
  ifndef SRCDIR
    SRCDIR=$(shell pwd)
  endif
  include $(TOPDIR)/src/bin/MakeVars
  vpath %.a $(libdir)
else
  SRCDIR=.
endif

CXXFLAGS += -I. -I$(SRCDIR)

RANLIB = ranlib
OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M

LIBCXXSRC = default_params.cc rr.cc dg.cc dgvertex.cc dgarc.cc gauss.cc \
oper.cc integral.cc iter.cc policy.cc strategy.cc policy_spec.cc flop.cc \
prefactors.cc context.cc memory.cc tactic.cc codeblock.cc dims.cc code.cc \
iface.cc class_registry.cc algebra.cc
LIBCXXOBJ = $(LIBCXXSRC:%.cc=%.$(OBJSUF))
LIBCXXDEP = $(LIBCXXSRC:%.cc=%.$(DEPSUF))
LIBOBJ = $(LIBCXXOBJ)
LIBTARGET = libINT.a

BUILDCXXSRC = build_libint.cc
BUILDCXXOBJ = $(BUILDCXXSRC:%.cc=%.$(OBJSUF))
BUILDCXXDEP = $(BUILDCXXSRC:%.cc=%.$(DEPSUF))

TESTCXXSRC = test.cc
TESTCXXOBJ = $(TESTCXXSRC:%.cc=%.$(OBJSUF))
TESTCXXDEP = $(TESTCXXSRC:%.cc=%.$(DEPSUF))

default:: build_libint

build_libint: $(BUILDCXXOBJ) $(LIBTARGET)
	$(CXX) -o $@ $(LDFLAGS) $^

depend:: $(TESTCXXDEP) $(BUILDCXXDEP) $(LIBCXXDEP)

test: $(TESTCXXOBJ) $(LIBTARGET) buildtest.o
	$(CXX) -o $@ $(LDFLAGS) $^

$(LIBTARGET): $(LIBOBJ)
	/bin/rm -f $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJ)
	$(RANLIB) $@

%.eps: %.ps
	ps2epsi $^ $@

%.pdf: %.ps
	ps2pdf13 $^

%.ps: %.dot
	dot -Tps -o $@ $^

graph.dot: test
	./$^

install::

clean:: oclean targetclean dependclean

oclean::
	-rm -f *.$(OBJSUF)

targetclean::
	-rm -f $(LIBTARGET) test build_libint

dependclean::
	-rm -f *.$(DEPSUF)

realclean:: clean

ifneq ($(CXXDEPENDSUF),none)
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

-include $(LIBCXXDEP)
-include $(BUILDCXXDEP)
-include $(TESTCXXDEP)
        
