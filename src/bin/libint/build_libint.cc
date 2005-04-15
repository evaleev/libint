
/**
  This program produces optimized source code to compute matrix elements (integrals)
  over Gaussian functions. Integrals are computed using recursive schemes. Generated source code
  is low-level C++ and can be used from other languages.
  */


/**
  The goals:
  1) Easier implementation of new RR :
     o) Number of centers != 4
     o) one-electron integrals
  2) Fortran and other interfaces from
  3) Vectorization
  4) General contractions
  */

//#include <libint_config.h>
#define ERI_MAX_AM 1
#define UNROLL_THRESH 1

#include <iostream>
#include <fstream>
#include <vector>
#include <rr.h>
#include <dg.h>
#include <typelist.h>
#include <integral.h>
#include <iter.h>
#include <policy_spec.h>
#include <intset_to_ints.h>
#include <strategy.h>

using namespace std;
using namespace libint2;

static int try_main (int argc, char* argv[]);

int main(int argc, char* argv[])
{
  try {
    try_main(argc,argv);
  }
  catch(std::exception& a) {
    cout << endl
         << "  WARNING! Caught a standard exception:" << endl
         << "    " << a.what() << endl << endl;
  }
}

static void print_header(std::ostream& os);
static void print_params(std::ostream& os);

static void build_TwoPRep_2b_2k(std::ostream& os,int lmax, int unroll_thresh);

int try_main (int argc, char* argv[])
{
  std::ostream& os = cout;
  
  print_header(os);
  print_params(os);
  
  int lmax_eri = ERI_MAX_AM;
  int unroll_thresh_eri = UNROLL_THRESH;
  
  build_TwoPRep_2b_2k(os,lmax_eri,unroll_thresh_eri);
}

void
print_header(std::ostream& os)
{
  os << "----------------------------------------------" << endl;
  os << " build_libint2: optimizing integrals compiler " << endl;
  os << "                by Edward F. Valeev           " << endl;
  os << "                and ideas by countless others " << endl;
  os << "----------------------------------------------" << endl << endl;
}

void
print_params(std::ostream& os)
{
  os << "ERI_MAX_AM    = " << static_cast<int>(ERI_MAX_AM) << endl;
  os << "UNROLL_THRESH = " << static_cast<int>(UNROLL_THRESH) << endl;
  os << endl;
}

void
build_TwoPRep_2b_2k(std::ostream& os, int lmax, int unroll_thresh)
{
  typedef TwoPRep_11_11<CGShell> TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {
          
          if (la+lb+lc+ld == 0)
            continue;
          
          SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
          SafePtr<Strategy> strat(new Strategy(unroll_thresh));
          SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
          os << "building " << abcd->description() << endl;
          SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(abcd);
          dg_xxxx->append_target(abcd_ptr);
          dg_xxxx->apply(strat);
          dg_xxxx->optimize_rr_out();
          dg_xxxx->traverse();
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
          
          SafePtr<CodeContext> context(new CppCodeContext());
          SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
          std::string prefix("libint_eri/");
          std::string decl_filename(prefix + context->label_to_name(abcd->label()));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(abcd->label()));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,abcd->label(),declfile,srcfile);
          
          os << "Max memory used = " << memman->max_memory_used() << endl;
          dg_xxxx->reset();
        }
      }
    }
  }
  
}

