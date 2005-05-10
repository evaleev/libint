
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

#include <libint_config.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <default_params.h>
#include <rr.h>
#include <dg.h>
#include <typelist.h>
#include <integral.h>
#include <iter.h>
#include <policy_spec.h>
#include <intset_to_ints.h>
#include <strategy.h>
#include <iface.h>
#include <r12kg12_11_11.h>
#include <vrr_11_r12kg12_11.h>

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
static void build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
static void build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);

int try_main (int argc, char* argv[])
{
  std::ostream& os = cout;
  
  // use default parameters
  SafePtr<CompilationParameters> cparams(new CompilationParameters);
  cparams->max_am_eri(1);
  // initialize code context to produce library API
  SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
  // make a list of computation labels
  Libint2Iface::Comps comps;
  comps.push_back("eri");
  comps.push_back("r12kg12");
  // intialize object to generate interface
  SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext,comps));
  
  print_header(os);
  cparams->print(os);
  
  build_TwoPRep_2b_2k(os,cparams,iface);
  build_R12kG12_2b_2k(os,cparams,iface);
  
  os << "Compilation finished. Goodbye." << endl;
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
build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  typedef TwoPRep_11_11<CGShell> TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am_eri();
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);
  
  iface->to_params(iface->define("MAX_AM_ERI",lmax));
  
  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
  SafePtr<Strategy> strat(new Strategy(cparams->unroll_threshold()));
  SafePtr<Tactic> tactic(new FirstChoiceTactic());
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {
          
          if (la+lb+lc+ld == 0)
            continue;
          
          SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
          os << "building " << abcd->description() << endl;
          SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(abcd);
          dg_xxxx->append_target(abcd_ptr);
          dg_xxxx->apply(strat,tactic);
          dg_xxxx->optimize_rr_out();
          dg_xxxx->traverse();
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
          
          SafePtr<CodeContext> context(new CppCodeContext(cparams));
          SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
          std::string prefix(cparams->source_directory());
          std::string decl_filename(prefix + context->label_to_name(abcd->label()));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(abcd->label()));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),abcd->label(),declfile,srcfile);
          
          ostringstream oss;
          oss << "  libint2_build_eri[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(abcd->label()))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());
          
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());
          
          os << "Max memory used = " << memman->max_memory_used() << endl;
          dg_xxxx->reset();
          declfile.close();
          srcfile.close();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
  
  //
  // generate explicit code for all recurrence relation that were not inlined
  //
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  std::string prefix(cparams->source_directory());
  dg_xxxx->generate_rr_code(context,prefix);
}


void
build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am_eri();
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);
  
  iface->to_params(iface->define("MAX_AM_R12kG12",lmax));
  
  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
  SafePtr<Strategy> strat(new Strategy(cparams->unroll_threshold()));
  SafePtr<Tactic> tactic(new FirstChoiceTactic());
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {
          
          if (la+lb+lc+ld == 0)
            continue;
          
          // k=0
          {
            typedef R12kG12_11_11<CGShell,0> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          // k=-1
          {
            typedef R12kG12_11_11<CGShell,-1> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          dg_xxxx->apply(strat,tactic);
          dg_xxxx->optimize_rr_out();
          dg_xxxx->traverse();
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
          
          std::string label;
          {
            ostringstream os;
            os << "(" << shells[la]->label() << " "
            << shells[lb]->label()
            << " | r_{12}^K * exp(-g*r_{12}^2) | "
            << shells[lc]->label() << " "
            << shells[ld]->label() << ")";
            label = os.str();
          }
          
          SafePtr<CodeContext> context(new CppCodeContext(cparams));
          SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
          std::string prefix(cparams->source_directory());
          std::string decl_filename(prefix + context->label_to_name(label));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(label));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),label,declfile,srcfile);
          
          ostringstream oss;
          oss << "  libint2_build_r12kg12[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());
          
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());
          
          os << "Max memory used = " << memman->max_memory_used() << endl;
          dg_xxxx->reset();
          declfile.close();
          srcfile.close();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
  
  //
  // generate explicit code for all recurrence relation that were not inlined
  //
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  std::string prefix(cparams->source_directory());
  dg_xxxx->generate_rr_code(context,prefix);
}

