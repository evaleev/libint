
#BOOSTPATH = /home/tool/evaleev/src/lib/boost
BOOSTPATH = /Users/evaleev/Development/General/src/lib/boost_1_32_0
CXX = gfilt
CXXFLAGS += -I. -g -I$(BOOSTPATH)

RANLIB = ranlib
OBJSUF = o

LIBCXXSRC = default_params.cc rr.cc dg.cc dgvertex.cc dgarc.cc gauss.cc \
oper.cc integral.cc iter.cc policy.cc strategy.cc policy_spec.cc flop.cc \
prefactors.cc context.cc memory.cc tactic.cc
LIBCXXOBJ = $(LIBCXXSRC:%.cc=%.$(OBJSUF))
LIBOBJ = $(LIBCXXOBJ)
LIBTARGET = libINT.a

BUILDCXXSRC = build_libint.cc
BUILDCXXOBJ = $(BUILDCXXSRC:%.cc=%.$(OBJSUF))

TESTCXXSRC = test.cc
TESTCXXOBJ = $(TESTCXXSRC:%.cc=%.$(OBJSUF))

default:: build_libint

build_libint: $(BUILDCXXOBJ) $(LIBTARGET)
	$(CXX) -o $@ $^

test: $(TESTCXXOBJ) $(LIBTARGET)
	$(CXX) -o $@ $^

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

clean:: oclean targetclean

oclean::
	-rm -f *.$(OBJSUF)

targetclean::
	-rm -f $(LIBTARGET) test build_libint
        
