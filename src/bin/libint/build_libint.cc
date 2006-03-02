
/**
  This program produces optimized source code to compute matrix elements (integrals)
  over Gaussian functions. Integrals are computed using recursive schemes. Generated source code
  is low-level C++ and can be used from other languages.

  Edward Valeev
  Atlanta/Oak Ridge
  December 2004.
  */

#define DO_TEST_ONLY 0

#include <libint2_config.h>

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
#include <vrr_11_r12kg12_11.h>
#include <r12kg12_11_11.h>
#include <tig12_11_11.h>
#include <graph_registry.h>

using namespace std;
using namespace libint2;

static void try_main (int argc, char* argv[]);

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
static void print_config(std::ostream& os);
// Put all configuration-specific API elements in here
static void config_to_api(const SafePtr<CompilationParameters>& cparams, SafePtr<Libint2Iface>& iface);
static void build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
static void build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
static void generate_rr_code(std::ostream& os, const SafePtr<CompilationParameters>& cparams);
static void test(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                 SafePtr<Libint2Iface>& iface);

void try_main (int argc, char* argv[])
{
  std::ostream& os = cout;
  
  // use default parameters
  SafePtr<CompilationParameters> cparams(new CompilationParameters);
  
#ifdef INCLUDE_ERI
  cparams->max_am_eri(ERI_MAX_AM);
  cparams->max_am_eri_opt(ERI_OPT_AM);
#endif
#ifdef INCLUDE_G12
  cparams->max_am_g12(G12_MAX_AM);
  cparams->max_am_g12_opt(G12_OPT_AM);
#endif
#if LIBINT_ENABLE_UNROLLING
  cparams->unroll_threshold(1000000000);
#endif
#ifdef LIBINT_VECTOR_LENGTH
  cparams->max_vector_length(LIBINT_VECTOR_LENGTH);
#endif
#ifdef LIBINT_VECTOR_METHOD
  {
    const std::string token("line");
    const bool vectorize_by_line = (token == LIBINT_VECTOR_METHOD);
    cparams->vectorize_by_line(vectorize_by_line);
  }
#endif
#if LIBINT_FLOP_COUNT
  cparams->count_flops(true);
#endif
#ifdef LIBINT_USER_DEFINED_FLOAT
  {
    const std::string realtype(LIBINT_USER_DEFINED_FLOAT);
    cparams->realtype(realtype);
  }
#endif
  
  // initialize code context to produce library API
  SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
  // make a list of computation labels
  Libint2Iface::Comps comps;
#ifdef INCLUDE_ERI
  comps.push_back("eri");
#endif
#ifdef INCLUDE_G12
  comps.push_back("r12kg12");
#endif
  // initialize object to generate interface
  SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext,comps));
  
  print_header(os);
  print_config(os);
  cparams->print(os);

#if !DO_TEST_ONLY
#ifdef INCLUDE_ERI
  build_TwoPRep_2b_2k(os,cparams,iface);
#endif
#ifdef INCLUDE_G12
  build_R12kG12_2b_2k(os,cparams,iface);
#endif
#endif

#if DO_TEST_ONLY
  test(os,cparams,iface);
#endif
  
  // Generate code for the set-level RRs
  generate_rr_code(os,cparams);
  
  // transfer some library configuration to library API
  config_to_api(cparams,iface);
  
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
print_config(std::ostream& os)
{
#ifdef INCLUDE_ERI
  os << "Will support ERI" << endl;
#endif
#ifdef INCLUDE_G12
  os << "Will support linear G12" << endl;
#endif
}

void
generate_rr_code(std::ostream& os, const SafePtr<CompilationParameters>& cparams)
{
  //
  // generate explicit code for all recurrence relation that were not inlined
  //
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  ImplicitDimensions::set_default_dims(cparams);
  std::string prefix(cparams->source_directory());

  SafePtr<RRStack> rrstack = RRStack::Instance();
  RRStack::citer_type it = rrstack->begin();
  while (it != rrstack->end()) {
    SafePtr<RecurrenceRelation> rr = (*it).second.second;
    std::string rrlabel = rr->label();
    os << "generating code for " << context->label_to_name(rrlabel) << " target=" << rr->rr_target()->label() << endl;
    
    std::string decl_filename(prefix + context->label_to_name(rrlabel));  decl_filename += ".h";
    std::string def_filename(prefix + context->label_to_name(rrlabel));  def_filename += ".cc";
    std::basic_ofstream<char> declfile(decl_filename.c_str());
    std::basic_ofstream<char> deffile(def_filename.c_str());
    
    rr->generate_code(context,ImplicitDimensions::default_dims(),declfile,deffile);
    
    declfile.close();
    deffile.close();
    
    // Remove RR to save resources
    rrstack->remove(rr);
    // next RR
    it = rrstack->begin();
  }
}

#ifdef INCLUDE_ERI
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
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {
          
          if (la+lb+lc+ld == 0)
            continue;
          if (la < lb || lc < ld || la+lb > lc+ld)
	    continue;
          
          // unroll only if max_am <= cparams->max_am_eri_opt()
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_eri_opt());
          dg_xxxx->registry()->can_unroll(need_to_optimize);
          dg_xxxx->registry()->do_cse(need_to_optimize);
          
          SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
          os << "building " << abcd->description() << endl;
          SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(abcd);
          dg_xxxx->append_target(abcd_ptr);
          dg_xxxx->apply(strat,tactic);
          dg_xxxx->optimize_rr_out();
          dg_xxxx->traverse();
#if DEBUG
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
#endif

          SafePtr<CodeContext> context(new CppCodeContext(cparams));
          SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
          std::string prefix(cparams->source_directory());
          std::string decl_filename(prefix + context->label_to_name(abcd->label()));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(abcd->label()));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),abcd->label(),declfile,srcfile);
          
          // update max stack size
          LibraryParameters& lparams = LibraryParameters::get_library_params();
          lparams.max_stack_size(memman->max_memory_used());
          
          ostringstream oss;
          oss << "  libint2_build_eri[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(abcd->label()))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());
          
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

#if DEBUG          
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxxx->reset();
          declfile.close();
          srcfile.close();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop  
}
#endif // INCLUDE_ERI

#ifdef INCLUDE_G12
void
build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am_g12();
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
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {

          if (la < lb || lc < ld || la+lb > lc+ld)
	    continue;
          bool ssss = false;
          if (la+lb+lc+ld == 0)
            ssss = true;
          
          // unroll only if max_am <= cparams->max_am_g12_opt()
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_g12_opt());
          dg_xxxx->registry()->can_unroll(need_to_optimize);
          dg_xxxx->registry()->do_cse(need_to_optimize);
          
          // k=0
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,0> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          // k=-1
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,-1> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          // [T_1,G12]
          if (true) {
            typedef TiG12_11_11<CGShell,0> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld]);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // [T_2,G12]
          if (true) {
            typedef TiG12_11_11<CGShell,1> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld]);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // k=2
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,2> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          dg_xxxx->apply(strat,tactic);
          dg_xxxx->optimize_rr_out();
          dg_xxxx->traverse();
#if DEBUG
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
#endif

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
          
          // update max stack size
          LibraryParameters& lparams = LibraryParameters::get_library_params();
          lparams.max_stack_size(memman->max_memory_used());

          ostringstream oss;
          oss << "  libint2_build_r12kg12[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());
          
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxxx->reset();
          declfile.close();
          srcfile.close();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop  
}
#endif // INCLUDE_G12

void
test(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
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
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  /*for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {*/
          {{{{
            int la=2;
            int lb=2;
            int lc=2;
            int ld=2;

          bool ssss = false;
          if (la+lb+lc+ld == 0)
            ssss = true;
          
          // k=0
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,0> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          // k=-1
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,-1> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          
          // [T_1,G12]
          if (true) {
            typedef TiG12_11_11<CGShell,0> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld]);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // [T_2,G12]
          if (true) {
            typedef TiG12_11_11<CGShell,1> int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld]);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // k=2
          if (!ssss) {
            typedef R12kG12_11_11<CGShell,2> int_type;
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
          
          // update max stack size
          LibraryParameters& lparams = LibraryParameters::get_library_params();
          lparams.max_stack_size(memman->max_memory_used());
          
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
}

void
config_to_api(const SafePtr<CompilationParameters>& cparams, SafePtr<Libint2Iface>& iface)
{
#ifdef INCLUDE_ERI
  iface->to_params(iface->define("SUPPORT_ERI",1));
  iface->to_params(iface->define("DERIV_ERI_ORDER",INCLUDE_ERI));
#endif
#ifdef INCLUDE_G12
  iface->to_params(iface->define("SUPPORT_G12",1));
  iface->to_params(iface->define("DERIV_G12_ORDER",INCLUDE_G12));
#endif
}

