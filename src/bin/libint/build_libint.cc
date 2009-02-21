
/**
  This program produces optimized source code to compute matrix elements (integrals)
  over Gaussian functions. Integrals are computed using recursive schemes. Generated source code
  is low-level C++ and can be used from other languages.

  Edward Valeev

  Atlanta/Oak Ridge (December 2004 - August 2006)
  Blacksburg (August 2006 - present)
  */

#include <libint2_config.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <default_params.h>
#include <rr.h>
#include <dg.h>
#include <dg.templ.h>
#include <integral.h>
#include <iter.h>
#include <policy_spec.h>
#include <strategy.h>
#include <iface.h>
#include <graph_registry.h>
#include <task.h>
#include <extract.h>
#include <dims.h>
#include <singl_stack.timpl.h>
#include <purgeable.h>

#include <master_ints_list.h>

using namespace std;
using namespace libint2;

#define STUDY_MEMORY_USAGE 0
long living_count = 0;

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
static void build_GenG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
			       SafePtr<Libint2Iface>& iface);
static void build_G12DKH_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
static void generate_rr_code(std::ostream& os, const SafePtr<CompilationParameters>& cparams);

void try_main (int argc, char* argv[])
{
  std::ostream& os = cout;

  // First must declare the tasks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
#ifdef INCLUDE_ERI
  taskmgr.add("eri");
#endif
#ifdef INCLUDE_G12
  taskmgr.add("r12kg12");
#endif
#ifdef INCLUDE_GENG12
  taskmgr.add("geng12");
#endif
#ifdef INCLUDE_G12DKH
  taskmgr.add("g12dkh");
#endif

  // use default parameters
  SafePtr<CompilationParameters> cparams(new CompilationParameters);

#ifdef INCLUDE_ERI
  cparams->max_am("eri",ERI_MAX_AM);
  cparams->max_am_opt("eri",ERI_OPT_AM);
#endif
#ifdef INCLUDE_G12
  cparams->max_am("r12kg12",G12_MAX_AM);
  cparams->max_am_opt("r12kg12",G12_OPT_AM);
#endif
#ifdef INCLUDE_GENG12
  cparams->max_am("geng12",GENG12_MAX_AM);
  cparams->max_am_opt("geng12",GENG12_OPT_AM);
#endif
#ifdef INCLUDE_G12DKH
  cparams->max_am("g12dkh",G12DKH_MAX_AM);
  cparams->max_am_opt("g12dkh",G12DKH_OPT_AM);
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
#if LIBINT_ACCUM_INTS
  cparams->accumulate_targets(true);
#else
  cparams->accumulate_targets(false);
#endif
#ifdef LIBINT_USER_DEFINED_FLOAT
  {
    const std::string realtype(LIBINT_USER_DEFINED_FLOAT);
    cparams->realtype(realtype);
  }
#endif
#ifdef LIBINT_API_PREFIX
  {
    const std::string api_prefix(LIBINT_API_PREFIX);
    cparams->api_prefix(api_prefix);
  }
#endif
#ifdef LIBINT_SINGLE_EVALTYPE
  {
    cparams->single_evaltype(true);
  }
#else
  throw std::runtime_error("Cannot generate specialized evaluator types yet");
#endif

  // initialize code context to produce library API
  SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
  // initialize object to generate interface
  SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext));

  print_header(os);
  print_config(os);
  // transfer some configuration parameters to the generated library API
  iface->to_params(iface->macro_define("CGSHELL_ORDERING",LIBINT_CGSHELL_ORDERING));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_STANDARD",LIBINT_CGSHELL_ORDERING_STANDARD));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_INTV3",LIBINT_CGSHELL_ORDERING_INTV3));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_GAMESS",LIBINT_CGSHELL_ORDERING_GAMESS));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_ORCA",LIBINT_CGSHELL_ORDERING_ORCA));
  cparams->print(os);

#ifdef INCLUDE_ERI
  build_TwoPRep_2b_2k(os,cparams,iface);
#endif
#ifdef INCLUDE_G12
  build_R12kG12_2b_2k(os,cparams,iface);
#endif
#ifdef INCLUDE_GENG12
  build_GenG12_2b_2k(os,cparams,iface);
#endif
#ifdef INCLUDE_G12DKH
  build_G12DKH_2b_2k(os,cparams,iface);
#endif

  // Generate code for the set-level RRs
  generate_rr_code(os,cparams);

  // print out the external symbols found for each task
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  for(tciter t=taskmgr.first(); t!=tend; ++t) {
    const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
    typedef TaskExternSymbols::SymbolList SymbolList;
    const SymbolList& symbols = tsymbols->symbols();
#if DEBUG
    // print out the labels
    std::cout << "Recovered labels for task " << t->label() << std::endl;
    typedef SymbolList::const_iterator citer;
    citer end = symbols.end();
    for(citer s=symbols.begin(); s!=end; ++s)
      std::cout << *s << std::endl;
#endif
  }

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
  os << "Will support G12 (commutators = ";
# if SUPPORT_T1G12
  os << "yes";
# else
  os << "false";
# endif
  os << ")" << endl;
#endif
#ifdef INCLUDE_GENG12
  os << "Will support generalized G12" << endl;
#endif
#ifdef INCLUDE_G12DKH
  os << "Will support G12DKH" << endl;
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
    std::string rrlabel = cparams->api_prefix() + rr->label();
    os << "generating code for " << context->label_to_name(rrlabel) << " target=" << rr->rr_target()->label() << endl;

    std::string decl_filename(prefix + context->label_to_name(rrlabel));  decl_filename += ".h";
    std::string def_filename(prefix + context->label_to_name(rrlabel));  def_filename += ".cc";
    std::basic_ofstream<char> declfile(decl_filename.c_str());
    std::basic_ofstream<char> deffile(def_filename.c_str());

    rr->generate_code(context,ImplicitDimensions::default_dims(),rrlabel,declfile,deffile);

    declfile.close();
    deffile.close();

    // Remove RR to save resources
    rrstack->remove(rr);
    // purge SingletonStacks, to save resources
    PurgeableStacks::Instance()->purge();
    // next RR
    it = rrstack->begin();
  }
}

#ifdef INCLUDE_ERI
void
build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  const std::string task("eri");
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_ERI",lmax));

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
  SafePtr<Strategy> strat(new Strategy());
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  //SafePtr<Tactic> tactic(new RandomChoiceTactic());
  //SafePtr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  for(int la=0; la<=lmax; la++) {
    for(int lb=0; lb<=lmax; lb++) {
      for(int lc=0; lc<=lmax; lc++) {
        for(int ld=0; ld<=lmax; ld++) {

          if (la+lb+lc+ld == 0)
            continue;

#if !GENERATE_FOR_ORCA
          if (la < lb || lc < ld || la+lb > lc+ld)
            continue;
#else
          if (la > lb || lc > ld || la+lb < lc+ld)
            continue;
#endif

#if STUDY_MEMORY_USAGE
          const int lim = 1;
          if (! (la == lim && lb == lim && lc == lim && ld == lim) )
            continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 1;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
	  dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
	  // Need to accumulate integrals?
	  dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

          SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],mType(0u));
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
	  std::string label(cparams->api_prefix() + abcd->label());
          std::string decl_filename(prefix + context->label_to_name(label));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(label));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),label,declfile,srcfile);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(memman->max_memory_used());
	  tparams->max_ntarget(1);

	  // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_eri[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

	  // need to declare this function internally
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

	  // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
	  // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
	  if (la == lmax &&
	      lb == lmax &&
	      lc == lmax &&
	      ld == lmax) {

	    extract_symbols(dg_xxxx);

	  }

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
  const std::string task("r12kg12");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_R12kG12",lmax));
#if SUPPORT_T1G12
  iface->to_params(iface->macro_define("SUPPORT_T1G12",1));
#else
  iface->to_params(iface->macro_define("SUPPORT_T1G12",0));
#endif

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
  SafePtr<Strategy> strat(new Strategy);
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
#if !SUPPORT_T1G12
	  if (ssss) continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 0;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
	  dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
	  // Need to accumulate integrals?
	  dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

      typedef R12kG12_11_11_sq int_type;
      typedef R12kG12 oper_type;
          // k=0
          if (!ssss) {
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(0));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // k=-1
          if (!ssss) {
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(-1));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

#if SUPPORT_T1G12
          // [T_1,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(0));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // [T_2,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(1));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
#endif

          // k=2
          if (!ssss) {
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(2));
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
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(memman->max_memory_used());
	  tparams->max_ntarget(5);

          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_r12kg12[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

	  // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
	  // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
	  if (la == lmax &&
	      lb == lmax &&
	      lc == lmax &&
	      ld == lmax) {

	    extract_symbols(dg_xxxx);

	  }

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

#if 0
#ifdef INCLUDE_GENG12
void
build_GenG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  const std::string task("geng12");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_GENG12",lmax));

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
  SafePtr<Strategy> strat(new Strategy);
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

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 1;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
	  dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
	  // Need to accumulate integrals?
	  dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

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

          // r1.r1 G12
          if (true) {
            typedef R1dotR1G12_11_11_sq int_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld]);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // r1.r2 G12
          if (true) {
            typedef R1dotR2G12_11_11_sq int_type;
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
            << " | exp(-a*r_1^2-a*r_2^2-g*r_{12}^2) | "
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
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(memman->max_memory_used());
	  tparams->max_ntarget(5);

          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_geng12[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

	  // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
	  // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
	  if (la == lmax &&
	      lb == lmax &&
	      lc == lmax &&
	      ld == lmax) {

	    extract_symbols(dg_xxxx);

	  }

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

#endif // INCLUDE_GENG12
#endif

#ifdef INCLUDE_G12DKH
void
build_G12DKH_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  const std::string task("g12dkh");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_G12DKH",lmax));

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
  SafePtr<Strategy> strat(new Strategy);
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

#if STUDY_MEMORY_USAGE
          const int lim = 5;
          if (! (la == lim && lb == lim && lc == lim && ld == lim) )
            continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la,lb),max(lc,ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 0;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
          dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
          // Need to accumulate integrals?
          dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

          // k=0
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(0));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // k=2
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(2));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // k=4
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(4));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // (G12prime.Div1)^2
          if (true) {
            typedef DivG12prime_xTx_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(0));
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // (G12prime.Div2)^2
          if (true) {
            typedef DivG12prime_xTx_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper_type(1));
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
            << " | r_{12}^(0,2,4) * exp(-g*r_{12}^2) | "
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
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(memman->max_memory_used());
          tparams->max_ntarget(3);

          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_g12dkh[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

      // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
      // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
      if (la == lmax &&
          lb == lmax &&
          lc == lmax &&
          ld == lmax) {

        extract_symbols(dg_xxxx);

      }

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

#endif // INCLUDE_G12DKH

void
config_to_api(const SafePtr<CompilationParameters>& cparams, SafePtr<Libint2Iface>& iface)
{
#ifdef INCLUDE_ERI
  iface->to_params(iface->macro_define("SUPPORT_ERI",1));
  iface->to_params(iface->macro_define("DERIV_ERI_ORDER",INCLUDE_ERI));
#endif
#ifdef INCLUDE_G12
  iface->to_params(iface->macro_define("SUPPORT_G12",1));
  iface->to_params(iface->macro_define("DERIV_G12_ORDER",INCLUDE_G12));
#endif
#ifdef INCLUDE_GENG12
  iface->to_params(iface->macro_define("SUPPORT_GENG12",1));
  iface->to_params(iface->macro_define("DERIV_GENG12_ORDER",INCLUDE_GENG12));
#endif
#ifdef INCLUDE_G12DKH
  iface->to_params(iface->macro_define("SUPPORT_G12DKH",1));
  iface->to_params(iface->macro_define("DERIV_G12DKH_ORDER",INCLUDE_G12DKH));
#endif
}

