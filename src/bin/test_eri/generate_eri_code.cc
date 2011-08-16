
#include <iostream>
#include <libint/util.h>
#include <libint/buildtest.h>
#include <libint/bfset.h>
#include <libint/master_ints_list.h>

using namespace libint2;

long living_count = 0;

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

  unsigned int am[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  typedef TwoPRep_11_11_sq ERIQtet;

  int try_main (int argc, char* argv[])
  {
#if LIBINT_CONTRACTED_INTS
    CGShell::set_contracted_default_value(true);
#else
    CGShell::set_contracted_default_value(false);
#endif
    TesterCmdLine<4> cmdline(argc,argv);
    const std::vector<unsigned int>& l = cmdline.am();

    // set this to the order of deriv ERI
    // will pick up from command line eventually
    const unsigned int deriv_order = 0;
    DerivIndexIterator<4u> diter(deriv_order);

    //
    // generate all targets
    //
    std::vector< SafePtr<ERIQtet> > targets;
    // iterate over the distribution of all derivative quanta among 12 different directions
    bool last_deriv = false;
    do {
      CGShell a(l[0]);  for(unsigned int xyz=0; xyz<3; ++xyz) a.deriv().inc(xyz, diter.value(0+xyz));
      CGShell b(l[1]);  for(unsigned int xyz=0; xyz<3; ++xyz) b.deriv().inc(xyz, diter.value(3+xyz));
      CGShell c(l[2]);  for(unsigned int xyz=0; xyz<3; ++xyz) c.deriv().inc(xyz, diter.value(6+xyz));
      CGShell d(l[3]);  for(unsigned int xyz=0; xyz<3; ++xyz) d.deriv().inc(xyz, diter.value(9+xyz));

      const SafePtr<ERIQtet> abcd0 = ERIQtet::Instance(a,b,c,d, 0u);
      targets.push_back(abcd0);
      last_deriv = diter.last();
      if (!last_deriv) diter.next();
    } while (!last_deriv);

    BuildTest<ERIQtet,true>(targets,
			    cmdline.size_to_unroll(),
			    cmdline.veclen(),
			    cmdline.vectorize_by_line(),
			    cmdline.do_cse(),
			    "eri0",
			    std::cout);

    return 0;
  }
};

