
#CXXFLAGS += -I. -g -I/home/tool/evaleev/src/lib/boost
CXXFLAGS += -I. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

CXXSRC = $(shell echo `ls *.cc`)
CXXOBJ = $(CXXSRC:%.cc=%.o)

test_eri: $(CXXOBJ)
	$(CXX) -o test_eri $^

clean::
	-rm -rf test_eri $(CXXOBJ)

