
CXXFLAGS += -I. -g -I~/src/lib/boost

test: test.o rr.o dg.o gauss.o oper.o integral.o iter.o
	$(CXX) -o test $^

clean::
	-rm -rf test.o rr.o dg.o gauss.o oper.o integral.o iter.o

