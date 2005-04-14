
#CXXFLAGS += -I. -I.. -g -I/home/tool/evaleev/src/lib/boost
#CXXFLAGS += -I. -I.. -O3 -xW -tpp7 -march=pentium4 -I/home/tool/evaleev/src/lib/boost
CXXFLAGS += -I. -I.. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

CXXGENSRC = $(shell echo `ls __*.cc`)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
LIBINT = ../libINT.a

TEST = test_eri
AUXSRC = prep_libint2.cc eri.cc
CXXTESTSRC = $(TEST).cc $(AUXSRC)
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)

TIME = time_eri
CXXTIMESRC = $(TIME).cc $(AUXSRC)
CXXTIMEOBJ = $(CXXTIMESRC:%.cc=%.o)

$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

clean::
	-rm -rf $(TEST) $(TIME) $(CXXTESTOBJ) $(CXXTIMEOBJ) $(CXXGENOBJ)

