
CXXFLAGS += -I. -I.. -g -I/home/tool/evaleev/src/lib/boost
#CXXFLAGS += -I. -I.. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

CXXSRC = $(shell echo `ls *.cc`)
CXXOBJ = $(CXXSRC:%.cc=%.o)
LIBINT = ../libINT.a

test_eri: $(CXXOBJ) $(LIBINT)
	$(CXX) -o $@ $^

clean::
	-rm -rf test_eri $(CXXOBJ)

