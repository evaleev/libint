
SRCDIR = ..
-include ../MakeVars

CPPFLAGS = -I. -I../include -DPREP_LIBINT2_SKIP_BOOST 
LDFLAGS = -L../lib -L../lib/.libs -lint2

TEST = test_eri

check:: $(TEST)
	./$(TEST) 0 3 && ./$(TEST) 1 2 && ./$(TEST) 2 1

$(TEST): test_eri.o eri.o
	$(CXX) $(LDFLAGS) -o $@ $^

clean:: targetclean
	-rm -rf *.o *.d

targetclean::
	-rm -rf $(TEST)

distclean:: clean

