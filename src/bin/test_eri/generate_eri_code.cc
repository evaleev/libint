
#include <iostream>
#include <libint/buildtest.h>
#include <libint/bfset.h>
#include <libint/master_ints_list.h>

using namespace libint2;

namespace {
  int try_main (int argc, char* argv[]);
};

int main (int argc, char* argv[])
{
  try {
    try_main(argc,argv);
  }
  catch(std::exception& a) {
    std::cerr << std::endl
         << "  WARNING! Caught a standard exception:" << std::endl
         << "    " << a.what() << std::endl << std::endl;
  }
}

namespace {

  unsigned int am[][1] = { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
  typedef TwoPRep_11_11_sq ERIQtet;

  int try_main (int argc, char* argv[])
  {
    TesterCmdLine<4> cmdline(argc,argv);
    const std::vector<unsigned int>& l = cmdline.am();
    const SafePtr<ERIQtet> quartet = ERIQtet::Instance(CGShell(am[l[0]]),
						       CGShell(am[l[1]]),
						       CGShell(am[l[2]]),
						       CGShell(am[l[3]]),
						       0u);
    SafePtr<DirectedGraph> pg1(new DirectedGraph);
    SafePtr<DirectedGraph> pg2(new DirectedGraph);
    SafePtr<DirectedGraph> pg3(new DirectedGraph);
    DirectedGraph g4;
    BuildTest<ERIQtet,true>(quartet,
			    cmdline.size_to_unroll(),
			    cmdline.veclen(),
			    cmdline.vectorize_by_line(),
			    cmdline.do_cse(),
			    "eri",
			    std::cout);
    return 0;
  }
};

