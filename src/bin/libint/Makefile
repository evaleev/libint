
CXXFLAGS += -I. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

test: test.o rr.o dg.o gauss.o oper.o integral.o iter.o policy.o strategy.o policy_spec.o
	$(CXX) -o test $^

clean::
	-rm -rf test.o rr.o dg.o gauss.o oper.o integral.o iter.o policy.o strategy.o policy_spec.o

