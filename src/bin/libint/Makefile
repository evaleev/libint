
CXXFLAGS += -I. -g

test: test.o rr.o dg.o gauss.o oper.o integral.o
	$(CXX) -o test $^

clean::
	-rm -rf test.o rr.o dg.o gauss.o oper.o integral.o

