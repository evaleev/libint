
CXXFLAGS += -I. -I.. -g -I/home/tool/evaleev/src/lib/boost
#CXXFLAGS += -I. -I.. -O3 -xW -tpp7 -march=pentium4 -mcpu=pentium4 -I/home/tool/evaleev/src/lib/boost
#CXXFLAGS += -I. -I.. -g -I/Users/evaleev/Development/General/src/lib/boost_1_32_0

CXXGENSRC = $(shell echo `ls __*.cc`)
CXXGENOBJ = $(CXXGENSRC:%.cc=%.o)
LIBINT = ../libINT.a
AUXSRC = prep_libint2v.cc eri.cc

TIME = time_eri_vec
CXXTIMESRC = $(TIME).cc $(AUXSRC)
CXXTIMEOBJ = $(CXXTIMESRC:%.cc=%.o)

$(TEST): $(CXXTESTOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

$(TIME): $(CXXTIMEOBJ) $(CXXGENOBJ) $(LIBINT)
	$(CXX) -o $@ $^

clean::
	-rm -rf $(TEST) $(TIME) $(CXXTESTOBJ) $(CXXTIMEOBJ) $(CXXGENOBJ)

