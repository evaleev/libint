
#CXXFLAGS += -I. -g -I/home/tool/evaleev/src/lib/boost
CXXFLAGS += -I. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

CXXSRC = test.cc rr.cc dg.cc gauss.cc oper.cc integral.cc iter.cc policy.cc \
strategy.cc policy_spec.cc flop.cc prefactors.cc context.cc
CXXOBJ = $(CXXSRC:%.cc=%.o)

test: $(CXXOBJ)
	$(CXX) -o test $^

graph.pdf: graph.ps
	ps2pdf13 graph.ps

graph.ps: graph.dot
	dot -Tps2 -o graph.ps graph.dot

graph.dot: test
	./test

clean::
	-rm -rf $(CXXOBJ)

