
CXXFLAGS += -I. -g -I/home/tool/evaleev/src/lib/boost

test: test.o rr.o dg.o gauss.o oper.o integral.o iter.o policy.o strategy.o policy_spec.o flop.o
	$(CXX) -o test $^

clean::
	-rm -rf test.o rr.o dg.o gauss.o oper.o integral.o iter.o policy.o strategy.o policy_spec.o flop.o

