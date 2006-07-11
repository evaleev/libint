
#include <iostream>
#include <libint2_config.h>
#include <libint/buildtest.h>
#include <libint/bfset.h>

using namespace libint2;

namespace {
  int try_main (int argc, char* argv[]);
  class PrintUsageAndDie : public std::runtime_error {
  public:
    PrintUsageAndDie(const char* what) : std::runtime_error(what) {}
  };
};

int main (int argc, char* argv[])
{
  try {
    try_main(argc,argv);
  }
  catch(PrintUsageAndDie& a) {
    std::cerr << std::endl
         << "  WARNING! Caught a standard exception:" << std::endl
         << "    " << a.what() << std::endl << std::endl;
    std::cout << "Usage: generate_eri_code a b c d size_to_unroll <vector_length> <vector_method>" << std::endl
	 << "       a,b,c,d -- angular momenta of functions in (ab|cd)" << std::endl
	 << "       size_to_unroll -- size of the largest quartet to be unrolled" << std::endl
	 << "       vector_length  -- (optional) max vector length. Defaults to 1." << std::endl << std::endl
	 << "       vector_method  -- (optional) vectorization method. Valid choices are 0 (by-block) and 1 (by-line). Defaults to 0." << std::endl
	 << "       do_cse  -- (optional) do Common Subexpression Elimination? Valid choices are 0 (no) and 1 (yes). Defaults to 0." << std::endl << std::endl;
  }
  catch(std::exception& a) {
    std::cerr << std::endl
         << "  WARNING! Caught a standard exception:" << std::endl
         << "    " << a.what() << std::endl << std::endl;
  }
}

namespace {

  const unsigned int max_am = 10;

  unsigned int am[][1] = { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
  CGShell sh_s(am[0]);
  CGShell sh_p(am[1]);
  CGShell sh_d(am[2]);
  CGShell sh_f(am[3]);
  CGShell sh_g(am[4]);
  CGShell sh_h(am[5]);
  CGShell sh_i(am[6]);
  CGShell sh_k(am[7]);
  CGShell sh_l(am[8]);
  CGShell sh_m(am[9]);
  CGShell sh_n(am[10]);
  typedef TwoPRep_11_11<CGShell> ERIQtet;

  int try_main (int argc, char* argv[])
  {
    if (argc < 6 || argc > 9)
      throw PrintUsageAndDie("incorrect number of command-line arguments");
    unsigned int la = atoi(argv[1]);
    unsigned int lb = atoi(argv[2]);
    unsigned int lc = atoi(argv[3]);
    unsigned int ld = atoi(argv[4]);
    unsigned int size_to_unroll = atoi(argv[5]);
    unsigned int veclen = 1;
    if (argc >= 7)
      veclen = atoi(argv[6]);
    bool vec_by_line = false;
    if (argc >= 8) {
      int vecmeth = atoi(argv[7]);
      if (vecmeth == 1)
        vec_by_line = true;
    }
    bool do_cse = false;
    if (argc >= 9) {
      int docse = atoi(argv[8]);
      if (docse == 1)
        do_cse = true;
    }
    
    if (la >= max_am || lb >= max_am || lc >= max_am || ld >= max_am)
      throw PrintUsageAndDie("Maximum angular momentum exceeded");
    const SafePtr<ERIQtet> quartet = ERIQtet::Instance(CGShell(am[la]),
						       CGShell(am[lb]),
						       CGShell(am[lc]),
						       CGShell(am[ld]),
						       0);

    std::cout << "generating code to compute " << quartet->label() << std::endl;
    std::cout << "size of ERI quartet representation = " << sizeof(ERIQtet) << std::endl;
    std::cout << "size of braket representation = " << sizeof(VectorBraket<CGShell>) << std::endl;
    std::cout << "size of DGVertex representation = " << sizeof(DGVertex) << std::endl;
    std::cout << "size of SafePtr<DGVertex> = " << sizeof(SafePtr<DGVertex>) << std::endl;

    // initialize cparams
    SafePtr<CompilationParameters> cparams(new CompilationParameters);
    cparams->max_am_eri(max_am);
    cparams->max_vector_length(veclen);
    cparams->vectorize_by_line(vec_by_line);
    cparams->count_flops(true);
#if LIBINT_ACCUM_INTS
    cparams->accumulate_targets(true);
#else
    cparams->accumulate_targets(false);
#endif
#ifdef LIBINT_API_PREFIX
  {
    const std::string api_prefix(LIBINT_API_PREFIX);
    cparams->api_prefix(api_prefix);
  }
#endif

    if (do_cse) {
      cparams->max_am_opt(max_am);
    }
    else {
      cparams->max_am_opt(0);
    }

    // set default dims
    ImplicitDimensions::set_default_dims(cparams);

    //SafePtr<StdRandomizePolicy> rpolicy(new StdRandomizePolicy(0.35));
    SafePtr<StdRandomizePolicy> rpolicy(new StdRandomizePolicy(0.00));
    SafePtr<Tactic> tactic(new FirstChoiceTactic<StdRandomizePolicy>(rpolicy));
    const SafePtr<MemoryManager> memman(new WorstFitMemoryManager);
    const std::string complabel("eri");
    BuildTest<ERIQtet,true>(quartet,cparams,size_to_unroll,std::cout,tactic,memman,complabel);

    return 0;
  }
};

