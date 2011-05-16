
#include <iostream>
#include <libint/buildtest.h>
#include <libint/bfset.h>

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
  return 0;
}

namespace {

  unsigned int am[][1] = { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
  typedef TiG12_11_11<CGShell,0> IntegralSet;

  int try_main (int argc, char* argv[])
  {
    TesterCmdLine<4> cmdline(argc,argv);
    std::vector<unsigned int> l = cmdline.am();
    const SafePtr<IntegralSet> integral = IntegralSet::Instance(CGShell(am[l[0]]),
								CGShell(am[l[1]]),
								CGShell(am[l[2]]),
								CGShell(am[l[3]]));
    BuildTest<IntegralSet,true>(integral,
				cmdline.size_to_unroll(),
				cmdline.veclen(),
				cmdline.vectorize_by_line(),
				cmdline.do_cse(),
				"t0g12",
				std::cout);
    return 0;
  }
};

