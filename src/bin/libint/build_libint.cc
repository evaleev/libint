
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
#include <integral.h>
#include <iter.h>
#include <policy_spec.h>
#include <strategy.h>
#include <iface.h>
#include <graph_registry.h>
#include <task.h>
#include <extract.h>
#include <dims.h>
#include <purgeable.h>
#include <buildtest.h>

#include <master_ints_list.h>

using namespace std;
using namespace libint2;

enum ShellSetType {
  ShellSetType_Standard = LIBINT_SHELL_SET_STANDARD,
  ShellSetType_ORCA     = LIBINT_SHELL_SET_ORCA
};
template <ShellSetType ShSet> struct ShellQuartetSetPredicate {
  // return true if this set of angular momenta is included
  static bool value(int la, int lb, int lc, int ld);
};
template <> struct ShellQuartetSetPredicate<ShellSetType_Standard> {
  static bool value(int la, int lb, int lc, int ld) {
    return la >= lb && lc >= ld && la+lb <= lc+ld;
  }
};
template <> struct ShellQuartetSetPredicate<ShellSetType_ORCA> {
  static bool value(int la, int lb, int lc, int ld) {
    return la <= lb && lc <= ld && ( la < lc || (la == lc && lb <= ld));
  }
};
template <ShellSetType ShSet> struct ShellTripletSetPredicate {
  // return true if this set of angular momenta is included
  static bool value(int lb, int lc, int cd);
};
template <> struct ShellTripletSetPredicate<ShellSetType_Standard> {
  static bool value(int lb, int lc, int ld) {
    return lc >= ld;
  }
};
template <> struct ShellTripletSetPredicate<ShellSetType_ORCA> {
  static bool value(int lb, int lc, int ld) {
    return lc <= ld;
  }
};

#define STUDY_MEMORY_USAGE 0
long living_count = 0;

std::string task_label(const std::string& prefix,
                       unsigned int deriv_level) {
  std::stringstream oss;
  if (deriv_level == 0)
    return prefix;
  else {
    oss << prefix << deriv_level;
    return oss.str();
  }
}

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
  return 0;
}

static void print_header(std::ostream& os);
static void print_config(std::ostream& os);
// Put all configuration-specific API elements in here
static void config_to_api(const SafePtr<CompilationParameters>& cparams, SafePtr<Libint2Iface>& iface);

#ifdef INCLUDE_ONEBODY
  static void build_OnePSep_1b_1k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                                 SafePtr<Libint2Iface>& iface);
  static void build_OnePNonSep_1b_1k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                                    SafePtr<Libint2Iface>& iface);
#endif

#ifdef INCLUDE_ERI
#define USE_GENERIC_ERI_BUILD 1
# if !USE_GENERIC_ERI_BUILD
static void build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
# else
static void build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface, unsigned int deriv_level);
# endif
#endif

#ifdef INCLUDE_ERI3
static void build_TwoPRep_1b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface, unsigned int deriv_level);
#endif

#ifdef INCLUDE_ERI2
static void build_TwoPRep_1b_1k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface, unsigned int deriv_level);
#endif

#ifdef INCLUDE_G12
static void build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
static void build_R12kG12_2b_2k_separate(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                         SafePtr<Libint2Iface>& iface);
#endif

#ifdef INCLUDE_G12DKH
static void build_G12DKH_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                                SafePtr<Libint2Iface>& iface);
#endif

void try_main (int argc, char* argv[])
{
  std::ostream& os = cout;

  // First must declare the tasks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.add("default");
#ifdef INCLUDE_ONEBODY
  for(unsigned int d=0; d<=INCLUDE_ONEBODY; ++d) {
    taskmgr.add( task_label("onebody",d) );
  }
#endif
#ifdef INCLUDE_ERI
  for(unsigned int d=0; d<=INCLUDE_ERI; ++d) {
    taskmgr.add( task_label("eri",d) );
  }
#endif
#ifdef INCLUDE_ERI3
  for(unsigned int d=0; d<=INCLUDE_ERI3; ++d) {
    taskmgr.add( task_label("3eri",d) );
  }
#endif
#ifdef INCLUDE_ERI2
  for(unsigned int d=0; d<=INCLUDE_ERI2; ++d) {
    taskmgr.add( task_label("2eri",d) );
  }
#endif
#ifdef INCLUDE_G12
   taskmgr.add("r12kg12");
# if !LIBINT_USE_COMPOSITE_EVALUATORS
   taskmgr.add("r12_0_g12");
   taskmgr.add("r12_2_g12");
# endif
#endif
#ifdef INCLUDE_GENG12
  taskmgr.add("geng12");
#endif
#ifdef INCLUDE_G12DKH
  taskmgr.add("g12dkh");
#endif

  // use default parameters
  SafePtr<CompilationParameters> cparams(new CompilationParameters);

  cparams->max_am("default",LIBINT_MAX_AM);
  cparams->max_am_opt("default",LIBINT_OPT_AM);
  cparams->num_bf("default",4);
#ifdef INCLUDE_ONEBODY
  for(unsigned int d=0; d<=0; ++d) {
    cparams->max_am( task_label("onebody", d) ,ONEBODY_MAX_AM);
    cparams->max_am_opt( task_label("onebody", d) ,ONEBODY_OPT_AM);
  }
#endif
#ifdef INCLUDE_ERI
  for(unsigned int d=0; d<=0; ++d) {
    cparams->max_am( task_label("eri", d) ,ERI_MAX_AM);
    cparams->max_am_opt( task_label("eri", d) ,ERI_OPT_AM);
  }
#endif
#ifdef INCLUDE_ERI3
  for(unsigned int d=0; d<=INCLUDE_ERI3; ++d) {
    cparams->max_am( task_label("3eri", d) ,ERI3_MAX_AM);
    cparams->max_am_opt( task_label("3eri", d) ,ERI3_OPT_AM);
  }
  for(unsigned int d=0; d<=INCLUDE_ERI3; ++d) {
    cparams->num_bf( task_label("3eri", d) ,3);
  }
#endif
#ifdef INCLUDE_ERI2
  for(unsigned int d=0; d<=INCLUDE_ERI2; ++d) {
    cparams->max_am( task_label("2eri", d) ,ERI2_MAX_AM);
    cparams->max_am_opt( task_label("2eri", d) ,ERI2_OPT_AM);
  }
  for(unsigned int d=0; d<=INCLUDE_ERI2; ++d) {
    cparams->num_bf( task_label("2eri", d) ,2);
  }
#endif
#ifdef INCLUDE_G12
# ifndef G12_MAX_AM
#   define LIBINT_MAX_AM G12_MAX_AM
# endif
# ifndef G12_OPT_AM
#   define LIBINT_OPT_AM G12_OPT_AM
# endif
  cparams->max_am("r12kg12",G12_MAX_AM);
  cparams->max_am_opt("r12kg12",G12_OPT_AM);
# if !LIBINT_USE_COMPOSITE_EVALUATORS
    cparams->max_am("r12_0_g12",G12_MAX_AM);
    cparams->max_am_opt("r12_0_g12",G12_OPT_AM);
    cparams->max_am("r12_2_g12",G12_MAX_AM);
    cparams->max_am_opt("r12_2_g12",G12_OPT_AM);
# endif
#endif
#ifdef INCLUDE_G12DKH
# ifndef G12DKH_MAX_AM
#   define LIBINT_MAX_AM G12DKH_MAX_AM
# endif
# ifndef G12DKH_OPT_AM
#   define LIBINT_OPT_AM G12DKH_OPT_AM
# endif
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
#if LIBINT_CONTRACTED_INTS
  cparams->contracted_targets(true);
  CGShell::set_contracted_default_value(true);
#else
  cparams->contracted_targets(false);
  CGShell::set_contracted_default_value(false);
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
  iface->to_params(iface->macro_define("USE_COMPOSITE_EVALUATORS",LIBINT_USE_COMPOSITE_EVALUATORS));
  iface->to_params(iface->macro_define("CARTGAUSS_MAX_AM",LIBINT_CARTGAUSS_MAX_AM));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING",LIBINT_CGSHELL_ORDERING));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_STANDARD",LIBINT_CGSHELL_ORDERING_STANDARD));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_INTV3",LIBINT_CGSHELL_ORDERING_INTV3));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_GAMESS",LIBINT_CGSHELL_ORDERING_GAMESS));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_ORCA",LIBINT_CGSHELL_ORDERING_ORCA));
  iface->to_params(iface->macro_define("SHELLQUARTET_SET",LIBINT_SHELL_SET));
  iface->to_params(iface->macro_define("SHELLQUARTET_SET_STANDARD",LIBINT_SHELL_SET_STANDARD));
  iface->to_params(iface->macro_define("SHELLQUARTET_SET_ORCA",LIBINT_SHELL_SET_ORCA));
  cparams->print(os);

#ifdef INCLUDE_ONEBODY
  build_OnePSep_1b_1k(os,cparams,iface,d);
  //build_OnePNonSep_2b_2k(os,cparams,iface,d);
#endif
#ifdef INCLUDE_ERI
# if !USE_GENERIC_ERI_BUILD
  build_TwoPRep_2b_2k(os,cparams,iface);
# else
  for(unsigned int d=0; d<=INCLUDE_ERI; ++d) {
    build_TwoPRep_2b_2k(os,cparams,iface,d);
  }
# endif
#endif
#ifdef INCLUDE_ERI3
  for(unsigned int d=0; d<=INCLUDE_ERI3; ++d) {
    build_TwoPRep_1b_2k(os,cparams,iface,d);
  }
#endif
#ifdef INCLUDE_ERI2
  for(unsigned int d=0; d<=INCLUDE_ERI2; ++d) {
    build_TwoPRep_1b_1k(os,cparams,iface,d);
  }
#endif
#ifdef INCLUDE_G12
# if LIBINT_USE_COMPOSITE_EVALUATORS
   build_R12kG12_2b_2k(os,cparams,iface);
# else
   build_R12kG12_2b_2k_separate(os,cparams,iface);
# endif
#endif
#ifdef INCLUDE_G12DKH
  build_G12DKH_2b_2k(os,cparams,iface);
#endif

  // Generate code for the set-level RRs
  std::deque<std::string> decl_filenames, def_filenames;
  generate_rr_code(os,cparams, decl_filenames, def_filenames);

#if DEBUG
  // print out the external symbols found for each task
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  for(tciter t=taskmgr.first(); t!=tend; ++t) {
    const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
    typedef TaskExternSymbols::SymbolList SymbolList;
    const SymbolList& symbols = tsymbols->symbols();
    // print out the labels
    std::cout << "Recovered labels for task " << t->label() << std::endl;
    typedef SymbolList::const_iterator citer;
    citer end = symbols.end();
    for(citer s=symbols.begin(); s!=end; ++s)
      std::cout << *s << std::endl;
  }
#endif

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
#ifdef INCLUDE_ONEBODY
  os << "Will support one-body integrals";
  if (INCLUDE_ONEBODY > 0)
    os << "(deriv order = " << INCLUDE_ONEBODY << ")";
  os << endl;
#endif
#ifdef INCLUDE_ERI
  os << "Will support ERI";
  if (INCLUDE_ERI > 0)
    os << "(deriv order = " << INCLUDE_ERI << ")";
  os << endl;
#endif
#ifdef INCLUDE_ERI3
  os << "Will support 3-center ERI";
  if (INCLUDE_ERI3 > 0)
    os << "(deriv order = " << INCLUDE_ERI3 << ")";
  os << endl;
#endif
#ifdef INCLUDE_ERI2
  os << "Will support 2-center ERI";
  if (INCLUDE_ERI2 > 0)
    os << "(deriv order = " << INCLUDE_ERI2 << ")";
  os << endl;
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
#ifdef INCLUDE_G12DKH
  os << "Will support G12DKH" << endl;
#endif
}

#ifdef INCLUDE_ONEBODY
void
build_OnePSep_1b_1k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface, unsigned int deriv_level)
{
  const std::string task = task_label("onebody", deriv_level);
  const std::string task_uc = task_label("ONEBODY", deriv_level);
  typedef OnePSep_1_1_sq OnePSep_sh_1_1;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(unsigned int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define( std::string("MAX_AM_") + task_uc,lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  SafePtr<DirectedGraph> dg(new DirectedGraph);
  SafePtr<Strategy> strat(new Strategy());
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int la=0; la<=lmax; la++) {
    for(unsigned int lb=0; lb<=lmax; lb++) {

          // skip s|s integrals -- no need to involve LIBINT here
          if (deriv_level == 0 && la == 0 && lb == 0)
            continue;

#if STUDY_MEMORY_USAGE
          const int lim = 1;
          if (! (la == lim && lb == lim) )
            continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(la,lb);
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 1;
          dg->registry()->unroll_threshold(unroll_threshold);
          dg->registry()->do_cse(need_to_optimize);
          dg->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
          // Need to accumulate integrals?
          dg->registry()->accumulate_targets(cparams->accumulate_targets());

          ////////////
          // loop over unique derivative index combinations
          ////////////
          // skip 1 center -- all derivatives with respect to that center can be
          // recovered using translational invariance conditions
          // which center to skip? -> A = 0, B = 1
          const unsigned int center_to_skip = 1;
          DerivIndexIterator<1> diter(deriv_level);
          std::vector< SafePtr<OnePSep_sh_1_1> > targets;
          bool last_deriv = false;
          do {
            CGShell a(la);
            CGShell b(lb);

            unsigned int center = 0;
            for(unsigned int i=0; i<2; ++i) {
              if (i == center_to_skip)
                continue;
              for(unsigned int xyz=0; xyz<3; ++xyz) {
                if (i == 0) a.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 1) b.deriv().inc(xyz, diter.value(3 * center + xyz));
              }
              ++center;
            }

            SafePtr<OnePSep_sh_1_1> abcd = OnePSep_sh_1_1::Instance(a,b,mType(0u));
            targets.push_back(abcd);
            last_deriv = diter.last();
            if (!last_deriv) diter.next();
          } while (!last_deriv);
          // append all derivatives as targets to the graph
          for(std::vector< SafePtr<OnePSep_sh_1_1> >::const_iterator t=targets.begin();
              t != targets.end();
              ++t) {
            SafePtr<DGVertex> t_ptr = dynamic_pointer_cast<DGVertex,OnePSep_sh_1_1>(*t);
            dg->append_target(t_ptr);
          }

          // make label that characterizes this set of targets
          // use the label of the nondifferentiated integral as a base
          std::string ab_label;
          {
            CGShell a(la);
            CGShell b(lb);
            SafePtr<OnePSep_sh_1_1> ab = OnePSep_sh_1_1::Instance(a,b,mType(0u));
            ab_label = ab->label();
          }
          // + derivative level (if deriv_level > 0)
          std::string label;
          {
            label = cparams->api_prefix();
            if (deriv_level != 0) {
              std::ostringstream oss;
              oss << "deriv" << deriv_level;
              label += oss.str();
            }
            label += ab_label;
          }

          std::string prefix(cparams->source_directory());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(targets.size());
          //os << " Max memory used = " << memman->max_memory_used() << std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_" << task << "[" << la << "][" << lb << "] = "
              << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

          // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
          // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
          if (la == lmax &&
              lb == lmax) {
            extract_symbols(dg);
          }

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg->reset();
          memman->reset();

        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
}
#endif

#ifdef INCLUDE_ERI
#if !USE_GENERIC_ERI_BUILD
void
build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  const std::string task("eri");
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(unsigned int l=0; l<=lmax; l++) {
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
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int la=0; la<=lmax; la++) {
    for(unsigned int lb=0; lb<=lmax; lb++) {
      for(unsigned int lc=0; lc<=lmax; lc++) {
        for(unsigned int ld=0; ld<=lmax; ld++) {

          if (la+lb+lc+ld == 0)
            continue;

          if (!ShellQuartetSetPredicate<static_cast<ShellSetType>(LIBINT_SHELL_SET)>::value(la,lb,lc,ld))
            continue;

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

          std::string prefix(cparams->source_directory());
          std::string label(cparams->api_prefix() + abcd->label());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(1);
          //os << " Max memory used = " << memman->max_memory_used() << std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_eri[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

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
          memman->reset();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
}
#else  // USE_GENERIC_ERI_BUILD
void
build_TwoPRep_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface, unsigned int deriv_level)
{
  const std::string task = task_label("eri", deriv_level);
  const std::string task_uc = task_label("ERI", deriv_level);
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(unsigned int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define( std::string("MAX_AM_") + task_uc,lmax));

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
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int la=0; la<=lmax; la++) {
    for(unsigned int lb=0; lb<=lmax; lb++) {
      for(unsigned int lc=0; lc<=lmax; lc++) {
        for(unsigned int ld=0; ld<=lmax; ld++) {

          // skip ss|ss integrals -- no need to involve LIBINT here
          if (deriv_level == 0 && la == 0 && lb == 0 && lc == 0 && ld == 0)
            continue;

          if (!ShellQuartetSetPredicate<static_cast<ShellSetType>(LIBINT_SHELL_SET)>::value(la,lb,lc,ld))
            continue;

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

          ////////////
          // loop over unique derivative index combinations
          ////////////
          // skip 1 center -- all derivatives with respect to that center can be
          // recovered using translational invariance conditions
          // which center to skip? -> A = 0, B = 1, C = 2, D = 3
          const unsigned int center_to_skip = 2;
          DerivIndexIterator<3> diter(deriv_level);
          std::vector< SafePtr<TwoPRep_sh_11_11> > targets;
          bool last_deriv = false;
          do {
            CGShell a(la);
            CGShell b(lb);
            CGShell c(lc);
            CGShell d(ld);

            unsigned int center = 0;
            for(unsigned int i=0; i<4; ++i) {
              if (i == center_to_skip)
                continue;
              for(unsigned int xyz=0; xyz<3; ++xyz) {
                if (i == 0) a.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 1) b.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 2) c.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 3) d.deriv().inc(xyz, diter.value(3 * center + xyz));
              }
              ++center;
            }

            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            targets.push_back(abcd);
            last_deriv = diter.last();
            if (!last_deriv) diter.next();
          } while (!last_deriv);
          // append all derivatives as targets to the graph
          for(std::vector< SafePtr<TwoPRep_sh_11_11> >::const_iterator t=targets.begin();
              t != targets.end();
              ++t) {
            SafePtr<DGVertex> t_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(*t);
            dg_xxxx->append_target(t_ptr);
          }

          // make label that characterizes this set of targets
          // use the label of the nondifferentiated integral as a base
          std::string abcd_label;
          {
            CGShell a(la);
            CGShell b(lb);
            CGShell c(lc);
            CGShell d(ld);
            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            abcd_label = abcd->label();
          }
          // + derivative level (if deriv_level > 0)
          std::string label;
          {
            label = cparams->api_prefix();
            if (deriv_level != 0) {
              std::ostringstream oss;
              oss << "deriv" << deriv_level;
              label += oss.str();
            }
            label += abcd_label;
          }

          std::string prefix(cparams->source_directory());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(targets.size());
          //os << " Max memory used = " << memman->max_memory_used() << std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_" << task << "[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

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
          memman->reset();

        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
}
# endif // USE_GENERIC_ERI_BUILD

#endif // INCLUDE_ERI

#ifdef INCLUDE_ERI3

void
build_TwoPRep_1b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface, unsigned int deriv_level)
{
  const std::string task = task_label("3eri", deriv_level);
  const std::string task_uc = task_label("3ERI", deriv_level);
  const std::string task_default("default");
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  const unsigned int lmax = cparams->max_am(task);
  const unsigned int lmax_default = deriv_level > 0 ? cparams->max_am(task_default) : lmax;
  for(unsigned int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define( std::string("MAX_AM_") + task_uc,lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  SafePtr<DirectedGraph> dg_xxx(new DirectedGraph);
  SafePtr<Strategy> strat(new Strategy());
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int lbra=0; lbra<=lmax; lbra++) {
    for(unsigned int lc=0; lc<=lmax_default; lc++) {
      for(unsigned int ld=0; ld<=lmax_default; ld++) {

          // skip ss|s integrals -- no need to involve LIBINT here
          if (deriv_level == 0 && lbra == 0 && lc == 0 && ld == 0)
            continue;

          // eliminate some cases depending on the desired convention
          if (!ShellTripletSetPredicate<static_cast<ShellSetType>(LIBINT_SHELL_SET)>::value(lbra,lc,ld))
            continue;

          // I will use 4-center recurrence relations and integrals, and have one center carry an s function
          // unfortunately, depending on the direction in which the build goes it must be A(0) or B(1)
          const unsigned int dummy_center = (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 0 : 1;

#if STUDY_MEMORY_USAGE
          const int lim = 1;
          if (! (lbra == lim && lc == lim && ld == lim) )
            continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(lc,ld),lbra);
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 1;
          dg_xxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxx->registry()->do_cse(need_to_optimize);
          dg_xxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
          // Need to accumulate integrals?
          dg_xxx->registry()->accumulate_targets(cparams->accumulate_targets());

          ////////////
          // loop over unique derivative index combinations
          ////////////
          // skip 1 center -- all derivatives with respect to that center can be
          // recovered using translational invariance conditions
          // which center to skip? -> the non-dummy one in bra -> A = 0, B = 1
          const unsigned int center_to_skip = (dummy_center == 0) ? 1 : 0;
          DerivIndexIterator<2> diter(deriv_level);
          std::vector< SafePtr<TwoPRep_sh_11_11> > targets;
          bool last_deriv = false;
          do {
            CGShell a((dummy_center == 0) ? 0 : lbra);
            CGShell b((dummy_center == 1) ? 0 : lbra);
            CGShell c(lc);
            CGShell d(ld);

            unsigned int center = 0;
            for(unsigned int i=0; i<4; ++i) {
              if (i == center_to_skip || i == dummy_center)
                continue;
              for(unsigned int xyz=0; xyz<3; ++xyz) {
                if (i == 0) a.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 1) b.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 2) c.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 3) d.deriv().inc(xyz, diter.value(3 * center + xyz));
              }
              ++center;
            }

            // use 4-center integrals
            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            targets.push_back(abcd);
            last_deriv = diter.last();
            if (!last_deriv) diter.next();
          } while (!last_deriv);
          // append all derivatives as targets to the graph
          for(std::vector< SafePtr<TwoPRep_sh_11_11> >::const_iterator t=targets.begin();
              t != targets.end();
              ++t) {
            SafePtr<DGVertex> t_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(*t);
            dg_xxx->append_target(t_ptr);
          }

          // make label that characterizes this set of targets
          // use the label of the nondifferentiated integral as a base
          std::string abcd_label;
          {
            CGShell a((dummy_center == 0) ? 0 : lbra);
            CGShell b((dummy_center == 1) ? 0 : lbra);
            CGShell c(lc);
            CGShell d(ld);
            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            abcd_label = abcd->label();
          }
          // + derivative level (if deriv_level > 0)
          std::string label;
          {
            label = cparams->api_prefix();
            if (deriv_level != 0) {
              std::ostringstream oss;
              oss << "deriv" << deriv_level;
              label += oss.str();
            }
            label += "eri3";
            label += abcd_label;
          }

          std::string prefix(cparams->source_directory());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg_xxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(targets.size());
          //os << " Max memory used = " << memman->max_memory_used() << std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_" << task << "[" << lbra << "][" << lc << "][" << ld << "] = "
              << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

          // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
          // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
          if (lbra == lmax &&
              lc == lmax &&
              ld == lmax) {
            extract_symbols(dg_xxx);
          }

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxx->reset();
          memman->reset();

      } // end of d loop
    } // end of c loop
  } // end of bra loop
}
#endif // INCLUDE_ERI3

#ifdef INCLUDE_ERI2

void
build_TwoPRep_1b_1k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface, unsigned int deriv_level)
{
  const std::string task = task_label("2eri", deriv_level);
  const std::string task_uc = task_label("2ERI", deriv_level);
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(unsigned int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define( std::string("MAX_AM_") + task_uc,lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  SafePtr<DirectedGraph> dg_xxx(new DirectedGraph);
  SafePtr<Strategy> strat(new Strategy());
  SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int lbra=0; lbra<=lmax; lbra++) {
    for(unsigned int lket=0; lket<=lmax; lket++) {

          // skip s|s integrals -- no need to involve LIBINT here
          if (deriv_level == 0 && lbra == 0 && lket == 0)
            continue;

          // I will use 4-center recurrence relations and integrals, and have two centers carry an s function
          // unfortunately, depending on the direction in which the build goes it must be A(0) and C(2) or B(1) and D(3)
          const unsigned int dummy_center1 = (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 0 : 1;
          const unsigned int dummy_center2 = (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 2 : 3;

#if STUDY_MEMORY_USAGE
          const int lim = 1;
          if (! (lbra == lim && lket == lim) )
            continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(lbra,lket);
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 1;
          dg_xxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxx->registry()->do_cse(need_to_optimize);
          dg_xxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
          // Need to accumulate integrals?
          dg_xxx->registry()->accumulate_targets(cparams->accumulate_targets());

          ////////////
          // loop over unique derivative index combinations
          ////////////
          // skip 1 center -- all derivatives with respect to that center can be
          // recovered using translational invariance conditions
          // which center to skip? -> the non-dummy one in bra -> A = 0, B = 1
          const unsigned int center_to_skip = (dummy_center1 == 0) ? 1 : 0;
          DerivIndexIterator<1> diter(deriv_level);
          std::vector< SafePtr<TwoPRep_sh_11_11> > targets;
          bool last_deriv = false;
          do {
            CGShell a((dummy_center1 == 0) ? 0 : lbra);
            CGShell b((dummy_center1 == 1) ? 0 : lbra);
            CGShell c((dummy_center2 == 2) ? 0 : lket);
            CGShell d((dummy_center2 == 3) ? 0 : lket);

            unsigned int center = 0;
            for(unsigned int i=0; i<4; ++i) {
              if (i == center_to_skip || i == dummy_center1 || i == dummy_center2)
                continue;
              for(unsigned int xyz=0; xyz<3; ++xyz) {
                if (i == 0) a.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 1) b.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 2) c.deriv().inc(xyz, diter.value(3 * center + xyz));
                if (i == 3) d.deriv().inc(xyz, diter.value(3 * center + xyz));
              }
              ++center;
            }

            // use 4-center integrals
            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            targets.push_back(abcd);
            last_deriv = diter.last();
            if (!last_deriv) diter.next();
          } while (!last_deriv);
          // append all derivatives as targets to the graph
          for(std::vector< SafePtr<TwoPRep_sh_11_11> >::const_iterator t=targets.begin();
              t != targets.end();
              ++t) {
            SafePtr<DGVertex> t_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_sh_11_11>(*t);
            dg_xxx->append_target(t_ptr);
          }

          // make label that characterizes this set of targets
          // use the label of the nondifferentiated integral as a base
          std::string abcd_label;
          {
            CGShell a((dummy_center1 == 0) ? 0 : lbra);
            CGShell b((dummy_center1 == 1) ? 0 : lbra);
            CGShell c((dummy_center2 == 2) ? 0 : lket);
            CGShell d((dummy_center2 == 3) ? 0 : lket);
            SafePtr<TwoPRep_sh_11_11> abcd = TwoPRep_sh_11_11::Instance(a,b,c,d,mType(0u));
            abcd_label = abcd->label();
          }
          // + derivative level (if deriv_level > 0)
          std::string label;
          {
            label = cparams->api_prefix();
            if (deriv_level != 0) {
              std::ostringstream oss;
              oss << "deriv" << deriv_level;
              label += oss.str();
            }
            label += "eri2";
            label += abcd_label;
          }

          std::string prefix(cparams->source_directory());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg_xxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size and # of targets
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(targets.size());
          //os << " Max memory used = " << memman->max_memory_used() << std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_" << task << "[" << lbra << "][" << lket << "] = "
              << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

          // For the most expensive (i.e. presumably complete) graph extract all precomputed quantities -- these will be members of the evaluator structure
          // also extract all RRs -- need to keep track of these to figure out which external symbols appearing in RR code belong to this task also
          if (lbra == lmax &&
              lket == lmax) {
            extract_symbols(dg_xxx);
          }

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxx->reset();
          memman->reset();

    } // end of ket loop
  } // end of bra loop
}
#endif // INCLUDE_ERI2

#ifdef INCLUDE_G12
void
build_R12kG12_2b_2k(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                    SafePtr<Libint2Iface>& iface)
{
  const std::string task("r12kg12");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for(unsigned int l=0; l<=lmax; l++) {
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
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

  for(unsigned int la=0; la<=lmax; la++) {
    for(unsigned int lb=0; lb<=lmax; lb++) {
      for(unsigned int lc=0; lc<=lmax; lc++) {
        for(unsigned int ld=0; ld<=lmax; ld++) {

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
            oper_type oper(0);
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // k=-1
          if (!ssss) {
            oper_type oper(-1);
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

#if SUPPORT_T1G12
          // [T_1,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            oper_type oper(0);
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // [T_2,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            oper_type oper(1);
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
#endif

          // [G12,[T1,G12]]
          if (!ssss) {
            typedef G12TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            oper_type oper(0); // doesn't matter whether T1 or T2 here
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
            os << "building " << abcd->description() << endl;
            SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          std::string _label;
          {
            ostringstream os;
            os << "(" << shells[la]->label() << " "
            << shells[lb]->label()
            << " | r_{12}^K * exp(-g*r_{12}^2) | "
            << shells[lc]->label() << " "
            << shells[ld]->label() << ")";
            _label = os.str();
          }

          std::string prefix(cparams->source_directory());
          std::string label(cparams->api_prefix() + _label);
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially generate code for its prerequisites
          GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames,
                       prefix, label, false);

          // update max stack size
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(5);

          ostringstream oss;
          oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_r12kg12[" << la << "][" << lb << "][" << lc << "]["
              << ld <<"] = " << context->label_to_name(label_to_funcname(label))
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
              i != decl_filenames.end();
              ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

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
          memman->reset();
        } // end of d loop
      } // end of c loop
    } // end of b loop
  } // end of a loop
}

#endif // INCLUDE_G12

#ifdef INCLUDE_G12
void
build_R12kG12_2b_2k_separate(std::ostream& os, const SafePtr<CompilationParameters>& cparams,
                             SafePtr<Libint2Iface>& iface)
{
  // do not support this if the commutator integrals are needed
#if SUPPORT_T1G12
  assert(false);
#endif

  // Note that because r12_-1_g12 integrals are evaluated using the same RR as ERI, no need to generate their code at all
  const int ntasks = 2;
  const char* task_names[] = {"r12_0_g12", "g12_T1_g12"};
  const char* task_NAMES[] = {"R12_0_R12", "G12_T1_G12"};

  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am("r12kg12");
  for(unsigned int l=0; l<=lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  for(int task=0; task<ntasks; ++task) {

    const std::string task_name(task_names[task]);

    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    taskmgr.current(task_name);
    iface->to_params(iface->macro_define(std::string("MAX_AM_") + task_NAMES[task],lmax));
    iface->to_params(iface->macro_define("SUPPORT_T1G12",0));

    SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
    SafePtr<Strategy> strat(new Strategy);
    SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
    SafePtr<CodeContext> context(new CppCodeContext(cparams));
    SafePtr<MemoryManager> memman(new WorstFitMemoryManager());

    for(unsigned int la=0; la<=lmax; la++) {
      for(unsigned int lb=0; lb<=lmax; lb++) {
        for(unsigned int lc=0; lc<=lmax; lc++) {
          for(unsigned int ld=0; ld<=lmax; ld++) {

            if (la+lb+lc+ld == 0)
              continue;

            if (!ShellQuartetSetPredicate<static_cast<ShellSetType>(LIBINT_SHELL_SET)>::value(la,lb,lc,ld))
              continue;

            using std::max;
            const unsigned int max_am = max(max(la,lb),max(lc,ld));
            const bool need_to_optimize = (max_am <= cparams->max_am_opt("r12kg12"));
            const unsigned int unroll_threshold = need_to_optimize ? cparams->unroll_threshold() : 0;
            dg_xxxx->registry()->unroll_threshold(unroll_threshold);
            dg_xxxx->registry()->do_cse(need_to_optimize);
            dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
            // Need to accumulate integrals?
            dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

            std::string _label;
            // k=0
            if (task == 0) {
              typedef R12kG12_11_11_sq int_type;
              typedef R12kG12 oper_type;
              oper_type oper(0);
#if LIBINT_CONTRACTED_INTS
              oper.descr().contract();
#endif
              SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
              os << "building " << abcd->description() << endl;
              SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
              dg_xxxx->append_target(abcd_ptr);
              _label = abcd_ptr->label();
            }

            // [G12,[T1,G12]]
            if (task == 1) {
              typedef G12TiG12_11_11_sq int_type;
              typedef int_type::OperType oper_type;
              oper_type oper(0); // doesn't matter whether T1 or T2 here
#if LIBINT_CONTRACTED_INTS
              oper.descr().contract();
#endif
              SafePtr<int_type> abcd = int_type::Instance(*shells[la],*shells[lb],*shells[lc],*shells[ld],0u,oper);
              os << "building " << abcd->description() << endl;
              SafePtr<DGVertex> abcd_ptr = dynamic_pointer_cast<DGVertex,int_type>(abcd);
              dg_xxxx->append_target(abcd_ptr);
              _label = abcd_ptr->label();
            }

            std::string prefix(cparams->source_directory());
            std::string label(cparams->api_prefix() + _label);
            std::deque<std::string> decl_filenames;
            std::deque<std::string> def_filenames;

            // this will generate code for this targets, and potentially generate code for its prerequisites
            GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                         decl_filenames, def_filenames,
                         prefix, label, false);

            // update max stack size
            const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
            tparams->max_stack_size(max_am, memman->max_memory_used());
            tparams->max_ntarget(1);

            ostringstream oss;
            oss << context->label_to_name(cparams->api_prefix()) << "libint2_build_" << task_names[task]
                << "[" << la << "][" << lb << "][" << lc << "]["
                << ld <<"] = " << context->label_to_name(label_to_funcname(label))
                << context->end_of_stat() << endl;
            iface->to_static_init(oss.str());

            // need to declare this function internally
            for(std::deque<std::string>::const_iterator i=decl_filenames.begin();
                i != decl_filenames.end();
                ++i) {
              oss.str("");
              oss << "#include <" << *i << ">" << endl;
              iface->to_int_iface(oss.str());
            }

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
            memman->reset();
          } // end of d loop
        } // end of c loop
      } // end of b loop
    } // end of a loop
  } // end of task loop
}

#endif // INCLUDE_G12

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

          SafePtr<CodeContext> context(new CppCodeContext(cparams));
          SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
          dg_xxxx->apply(strat,tactic);
          dg_xxxx->optimize_rr_out(context);
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
          std::string prefix(cparams->source_directory());
          std::string decl_filename(prefix + context->label_to_name(label));  decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(label));  src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),label,declfile,srcfile);

          // update max stack size
          const SafePtr<TaskParameters>& tparams = taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
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
#ifdef INCLUDE_ERI3
  iface->to_params(iface->macro_define("SUPPORT_ERI3",1));
  iface->to_params(iface->macro_define("DERIV_ERI3_ORDER",INCLUDE_ERI3));
#endif
#ifdef INCLUDE_ERI2
  iface->to_params(iface->macro_define("SUPPORT_ERI2",1));
  iface->to_params(iface->macro_define("DERIV_ERI2_ORDER",INCLUDE_ERI2));
#endif
#ifdef INCLUDE_G12
  iface->to_params(iface->macro_define("SUPPORT_G12",1));
  iface->to_params(iface->macro_define("DERIV_G12_ORDER",INCLUDE_G12));
#endif
#ifdef INCLUDE_G12DKH
  iface->to_params(iface->macro_define("SUPPORT_G12DKH",1));
  iface->to_params(iface->macro_define("DERIV_G12DKH_ORDER",INCLUDE_G12DKH));
#endif
}

