
CXXFLAGS += -I. -g

test: test.o rr.o dg.o gauss.o oper.o integral.o iter.o traits.o
	$(CXX) -o test $^

clean::
	-rm -rf test.o rr.o dg.o gauss.o oper.o integral.o iter.o traits.o

