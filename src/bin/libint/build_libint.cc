/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
  This program produces optimized source code to compute matrix elements
  (integrals) over Gaussian functions. Integrals are computed using recursive
  schemes. Generated source code is low-level C++ and can be used from other
  languages.

  Edward Valeev

  Atlanta/Oak Ridge (December 2004 - August 2006)
  Blacksburg (August 2006 - present)
  */

#include <boost/preprocessor.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#if not BOOST_PP_VARIADICS  // no variadic macros? your compiler is out of date!
                            // (should not be possible since variadic macros are
                            // part of C++11)
#error \
    "your compiler does not provide variadic macros (but does support C++11), something is seriously broken, please create an issue at https://github.com/evaleev/libint/issues"
#endif

#include <buildtest.h>
#include <default_params.h>
#include <dg.h>
#include <dims.h>
#include <extract.h>
#include <global_macros.h>
#include <graph_registry.h>
#include <iface.h>
#include <integral.h>
#include <iter.h>
#include <libint2/cgshell_ordering.h>
#include <libint2/config.h>
#include <libint2/deriv_iter.h>
#include <libint2/shgshell_ordering.h>
#include <master_ints_list.h>
#include <policy_spec.h>
#include <purgeable.h>
#include <rr.h>
#include <strategy.h>
#include <task.h>

using namespace std;
using namespace libint2;

CodeGenProgress g_progress;

enum ShellSetType {
  ShellSetType_Standard = LIBINT_SHELL_SET_STANDARD,
  ShellSetType_ORCA = LIBINT_SHELL_SET_ORCA
};
template <ShellSetType ShSet>
struct ShellQuartetSetPredicate {
  // return true if this set of angular momenta is included
  static bool value(int la, int lb, int lc, int ld, bool p1p2_swappable = true);
};

/**
 * standard ordering for angular momenta la, lb, lc, ld
 * @param p1p2_swappable whether operator allows swaps of particle 1 and 2
 * functions (e.g., not allowed for Coulombσpσp but allowed
 * for Coulomb (TwoPRep)).
 * @param bra_ket_coswappable whether need to swap within both bra and ket.
 * Not individually swapping of either ket of bra allowed
 * ( e.g., for σpσpCoulombσpσp)
 */
template <>
struct ShellQuartetSetPredicate<ShellSetType_Standard> {
  static bool value(int la, int lb, int lc, int ld, bool p1p2_swappable = true,
                    bool bra_ket_coswappable = false) {
    if (bra_ket_coswappable)
      return lc >= ld && (la + lb < lc + ld ||
                          (la + lb == lc + ld && std::max(la, lb) <= lc));
    else
      return la >= lb && lc >= ld && (!p1p2_swappable || la + lb <= lc + ld);
  }
};
template <>
struct ShellQuartetSetPredicate<ShellSetType_ORCA> {
  static bool value(int la, int lb, int lc, int ld, bool p1p2_swappable = true,
                    bool bra_ket_coswappable = false) {
    if (bra_ket_coswappable)
      return lc <= ld && (la + lb > lc + ld ||
                          (la + lb == lc + ld && std::min(la, lb) >= lc));
    else
      return la <= lb && lc <= ld &&
             (!p1p2_swappable || (la < lc || (la == lc && lb <= ld)));
  }
};
template <ShellSetType ShSet>
struct ShellTripletSetPredicate {
  // return true if this set of angular momenta is included
  static bool value(int lb, int lc, int cd);
};
template <>
struct ShellTripletSetPredicate<ShellSetType_Standard> {
  static bool value(int lb, int lc, int ld) { return lc >= ld; }
};
template <>
struct ShellTripletSetPredicate<ShellSetType_ORCA> {
  static bool value(int lb, int lc, int ld) { return lc <= ld; }
};

#define STUDY_MEMORY_USAGE 0
long living_count = 0;

size_t l_to_cgshellsize(size_t l) { return (l + 1) * (l + 2) / 2; }

std::string task_label(const std::string& prefix, unsigned int deriv_level) {
  std::stringstream oss;
  if (deriv_level == 0)
    return prefix;
  else {
    oss << prefix << deriv_level;
    return oss.str();
  }
}

/**
 * Returns n-th token. If n > number of tokens, returns the last one.
 * @param c_str input string
 * @param delimiter delimiter/separator character
 * @param n
 * @return
 */
template <typename T>
T token(const char* c_str, char delimiter, std::size_t n = 0) {
  T result;
  std::string result_str;

  // replace all occurences of delimiter in str with whitespaces
  std::string str(c_str);
  std::cout << "token<>: str=" << str << std::endl;
  std::string::size_type pos;
  while ((pos = str.find(delimiter)) != std::string::npos) {
    str[pos] = ' ';
  }

  // tokenize
  std::istringstream iss(str);
  for (size_t i = 0; i <= n && !iss.eof(); ++i) iss >> result;

  return result;
}

/**
 * Returns the number of tokens.
 * @param c_str input string
 * @param delimiter delimiter/separator character
 * @return the number of tokens
 */
inline size_t ntokens(const char* c_str, char delimiter) {
  size_t n = 1;

  // replace all occurences of delimiter in str with whitespaces
  std::string str(c_str);
  std::cout << "ntokens: str=" << str << std::endl;
  std::string::size_type pos;
  while ((pos = str.find(delimiter)) != std::string::npos) {
    str[pos] = ' ';
    ++n;
  }

  return n;
}

static void try_main(int argc, char* argv[]);

int main(int argc, char* argv[]) {
  int return_code = 0;
  try {
    try_main(argc, argv);
  } catch (std::exception& a) {
    cout << endl
         << "  WARNING! Caught a standard exception:" << endl
         << "    " << a.what() << endl
         << endl;
    return_code = 1;
  } catch (...) {
    cout << endl << "  WARNING! Caught an unknown exception" << endl << endl;
    return_code = 1;
  }

  return return_code;
}

static void print_header(std::ostream& os);
static void print_config(std::ostream& os);
// Put all configuration-specific API elements in here
static void config_to_api(const std::shared_ptr<CompilationParameters>& cparams,
                          std::shared_ptr<Libint2Iface>& iface);

#ifdef LIBINT_INCLUDE_ERI
#define USE_GENERIC_ERI_BUILD 1
#endif
#if defined(LIBINT_INCLUDE_ERI) || defined(LIBINT_INCLUDE_RKB_ERI)
#if defined(USE_GENERIC_ERI_BUILD) && !USE_GENERIC_ERI_BUILD
template <typename OperType>
static void build_TwoPRep_2b_2k(
    std::ostream& os, std::string label,
    const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface);
#else
template <typename OperType>
static void build_TwoPRep_2b_2k(
    std::ostream& os, std::string label,
    const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface, unsigned int deriv_level = 0);
#endif
#endif

#ifdef LIBINT_INCLUDE_ERI3
static void build_TwoPRep_1b_2k(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface, unsigned int deriv_level);
#endif

#ifdef LIBINT_INCLUDE_ERI2
static void build_TwoPRep_1b_1k(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface, unsigned int deriv_level);
#endif

#ifdef LIBINT_INCLUDE_G12
static void build_R12kG12_2b_2k(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface);
static void build_R12kG12_2b_2k_separate(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface);
#endif

#ifdef LIBINT_INCLUDE_G12DKH
static void build_G12DKH_2b_2k(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface);
#endif

#ifdef LIBINT_INCLUDE_ONEBODY

#if LIBINT_SUPPORT_ONEBODYINTS == 0
#error \
    "change LIBINT_SUPPORT_ONEBODYINTS in global_macros.h to 1 if need 1-body ints"
#endif

namespace {
template <typename OperType>
struct AuxQuantaType;
template <>
struct AuxQuantaType<ElecPotOper> {
  typedef mType type;
};
template <typename OperType>
struct AuxQuantaType {
  typedef EmptySet type;
};

}  // namespace

template <typename OperDescrType>
OperDescrType make_descr(int, int = 0, int = 0) {
  return OperDescrType();
}
template <>
CartesianMultipole_Descr<3u> make_descr<CartesianMultipole_Descr<3u>>(int x,
                                                                      int y,
                                                                      int z) {
  CartesianMultipole_Descr<3u> result;
  // cartesian multipole quanta
  result.inc(0, x);
  result.inc(1, y);
  result.inc(2, z);
  return result;
}
template <>
SphericalMultipole_Descr make_descr<SphericalMultipole_Descr>(int l, int m,
                                                              int) {
  SphericalMultipole_Descr result(l, m);
  return result;
}
template <>
σpVσp_Descr make_descr<σpVσp_Descr>(int p, int, int) {
  return σpVσp_Descr(p);
}

template <>
Coulombσpσp_Descr make_descr<Coulombσpσp_Descr>(int p, int, int) {
  return Coulombσpσp_Descr(p);
}

template <>
σpσpCoulombσpσp_Descr make_descr<σpσpCoulombσpσp_Descr>(int p, int, int) {
  return σpσpCoulombσpσp_Descr(p);
}

template <typename _OperType>
void build_onebody_1b_1k(std::ostream& os, std::string label,
                         const std::shared_ptr<CompilationParameters>& cparams,
                         std::shared_ptr<Libint2Iface>& iface,
                         unsigned int deriv_level) {
  // implement overlap as a special case of cartesian emultipole
  using OperType =
      typename std::conditional<std::is_same<_OperType, OverlapOper>::value,
                                CartesianMultipoleOper<3u>, _OperType>::type;
  const std::string task = task_label(label, deriv_level);
  typedef CGShell BFType;
  typedef typename OperType::Descriptor OperDescrType;
  typedef GenIntegralSet_1_1<CGShell, OperType,
                             typename AuxQuantaType<OperType>::type>
      Onebody_sh_1_1;

  vector<BFType*> shells;
  unsigned int lmax = cparams->max_am(task);
  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new BFType(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define(std::string("MAX_AM_") + task, lmax));

  const auto nullaux = typename Onebody_sh_1_1::AuxIndexType(0u);

  // optionally skip derivative property ints
#ifdef LIBINT_DISABLE_ONEBODY_PROPERTY_DERIVS
  const auto property_operator =
      !(std::is_same<_OperType, OverlapOper>::value ||
        std::is_same<_OperType, KineticOper>::value ||
        std::is_same<_OperType, ElecPotOper>::value);
  if (property_operator && deriv_level > 0) return;
#endif
  // derivatives of spherical multipole integrals are not implemented
  {
    if (std::is_same<_OperType, SphericalMultipoleOper>::value &&
        deriv_level > 0)
      return;
    // throw std::invalid_argument(
    //     "derivatives of spherical multipole ints are not yet implemented");
  }

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  std::shared_ptr<DirectedGraph> dg(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy());
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

  for (unsigned int la = 0; la <= lmax; la++) {
    for (unsigned int lb = 0; lb <= lmax; lb++) {
      // skip s|s overlap and elecpot integrals -- no need to involve LIBINT
      // here
      if (deriv_level == 0 && la == 0 && lb == 0 &&
          (std::is_same<_OperType, OverlapOper>::value ||
           std::is_same<_OperType, ElecPotOper>::value))
        continue;

      std::shared_ptr<Tactic> tactic(new TwoCenter_OS_Tactic(la, lb));

      // this will hold all target shell sets
      std::vector<std::shared_ptr<Onebody_sh_1_1>> targets;

      /////////////////////////////////
      // loop over operator components
      /////////////////////////////////
      // most important operators have 1 component ...
      std::vector<OperDescrType> descrs(1);  // operator descriptors
      // important EXCEPTION: multipole moments
      if (std::is_same<_OperType, CartesianMultipoleOper<3u>>::value) {
        // reset descriptors array
        descrs.resize(0);

        // parse the label ... 1emultipole means include multipoles of order 0
        // (overlap) and 1 (dipole)
        unsigned int max_multipole_order = 0;
        assert(label != "overlap");
        auto key_pos = label.find("emultipole");
        assert(key_pos != std::string::npos);
        std::string tmp = label;
        tmp.erase(key_pos, std::string::npos);
        istringstream iss(tmp);
        iss >> max_multipole_order;
        assert(max_multipole_order > 0);

        // iterate over operators and construct their descriptors
        for (int multipole_order = 0; multipole_order <= max_multipole_order;
             ++multipole_order) {
          // we iterate over them same way as over cartesian Gaussian shells
          int x, y, z;
          FOR_CART(x, y, z, multipole_order)
          descrs.push_back(make_descr<OperDescrType>(x, y, z));
          END_FOR_CART
        }
      }
      if (std::is_same<_OperType, SphericalMultipoleOper>::value) {
        // reset descriptors array
        descrs.resize(0);
        // iterate over operators and construct their descriptors
        for (int l = 0; l <= LIBINT_MULTIPOLE_MAX_ORDER; ++l) {
          // we iterate over them using the *standard* solid harmonic Gaussian
          // ordering
          int m;
          FOR_SOLIDHARM_STANDARD(l, m)
          descrs.push_back(make_descr<OperDescrType>(l, m));
          END_FOR_SOLIDHARM
        }
      }
      if (std::is_same<_OperType, σpVσpOper>::value) {
        // reset descriptors array
        descrs.resize(0);
        // iterate over pauli components
        for (int p = 0; p != 4; ++p) {
          descrs.emplace_back(make_descr<OperDescrType>(p));
        }
      }

      // derivative index is the outermost (slowest running)
      // operator component is second slowest

      ////////////
      // loop over unique derivative index combinations
      ////////////
      // NB translational invariance is now handled by CR_DerivGauss
      CartesianDerivIterator<2> diter(deriv_level);
      bool last_deriv = false;
      do {
        BFType a(la);
        BFType b(lb);

        for (unsigned int c = 0; c != 2; ++c) {
          const unsigned int ndir =
              std::is_same<BFType, CGShell>::value ? 3 : 1;
          for (unsigned int xyz = 0; xyz < ndir; ++xyz) {
            if (c == 0) a.deriv().inc(xyz, (*diter).at(xyz));
            if (c == 1) b.deriv().inc(xyz, (*diter).at(3 + xyz));
          }
        }

        // operator component loop
        for (unsigned int op = 0; op != descrs.size(); ++op) {
          OperType oper(descrs[op]);

          std::shared_ptr<Onebody_sh_1_1> target =
              Onebody_sh_1_1::Instance(a, b, nullaux, oper);
          targets.push_back(target);
        }  // loop over operator components

        last_deriv = diter.last();
        if (!last_deriv) diter.next();
      } while (!last_deriv);  // loop over derivatives

      // unroll only if max_am <= cparams->max_am_opt(task)
      using std::max;
      const unsigned int max_am = max(la, lb);
      const auto nopers = descrs.size();
      const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
      // decide whether to unroll based on the aggregate size of
      // undifferentiated quartets
      const bool need_to_unroll =
          nopers * l_to_cgshellsize(la) * l_to_cgshellsize(lb) <=
          cparams->unroll_threshold();
      const unsigned int unroll_threshold =
          need_to_optimize && need_to_unroll
              ? std::numeric_limits<unsigned int>::max()
              : 0;

      dg->registry()->unroll_threshold(unroll_threshold);
      dg->registry()->do_cse(need_to_optimize);
      dg->registry()->condense_expr(condense_expr(
          cparams->unroll_threshold(), cparams->max_vector_length() > 1));
      // Need to accumulate integrals?
      dg->registry()->accumulate_targets(cparams->accumulate_targets());

      // shove all targets on the graph, IN ORDER
      for (auto t = targets.begin(); t != targets.end(); ++t) {
        std::shared_ptr<DGVertex> t_ptr =
            std::dynamic_pointer_cast<DGVertex, Onebody_sh_1_1>(*t);
        dg->append_target(t_ptr);
      }

      // make label that characterizes this set of targets
      std::string eval_label;
      {
        std::ostringstream oss;
        oss << "_" << label;
        if (deriv_level > 0) oss << "deriv" << deriv_level;
        BFType a(la);
        BFType b(lb);
        oss << "_" << a.label() << "_" << b.label();
        eval_label = oss.str();
      }

      g_progress.current_task = eval_label;
      g_progress.print();

      std::string prefix(cparams->source_directory());
      std::deque<std::string> decl_filenames;
      std::deque<std::string> def_filenames;

      // this will generate code for this targets, and potentially generate code
      // for its prerequisites
      GenerateCode(dg, context, cparams, strat, tactic, memman, decl_filenames,
                   def_filenames, prefix, eval_label, false);

      // update max stack size and # of targets
      const std::shared_ptr<TaskParameters>& tparams =
          taskmgr.current().params();
      tparams->max_stack_size(max_am, memman->max_memory_used());
      tparams->max_ntarget(targets.size());
      // os << " Max memory used = " << memman->max_memory_used() << std::endl;

      // set pointer to the top-level evaluator function
      ostringstream oss;
      oss << context->label_to_function_name("libint2_build_" + task) << "["
          << la << "][" << lb
          << "] = " << context->label_to_function_name(eval_label)
          << context->end_of_stat() << endl;
      iface->to_static_init(oss.str());

      // need to declare this function internally
      for (std::deque<std::string>::const_iterator i = decl_filenames.begin();
           i != decl_filenames.end(); ++i) {
        oss.str("");
        oss << "#include <" << *i << ">" << endl;
        iface->to_int_iface(oss.str());
      }

#if DEBUG
      os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
      dg->reset();
      memman->reset();

    }  // end of b loop
  }    // end of a loop
}
#endif

void try_main(int argc, char* argv[]) {
  std::ostream& os = cout;

  // First must declare the tasks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.add("default");
#if defined(LIBINT_MAX_AM_LIST)
  for (unsigned int d = 1; d < ntokens(LIBINT_MAX_AM_LIST, ','); ++d) {
    taskmgr.add(task_label("default", d));
  }
#endif
#ifdef LIBINT_INCLUDE_ONEBODY

// overlap, kinetic, elecpot cannot be omitted
#define BOOST_PP_ONEBODY_TASK_TUPLE                                  \
  (overlap, kinetic, elecpot, 1emultipole, 2emultipole, 3emultipole, \
   sphemultipole, opVop)
#define BOOST_PP_ONEBODY_TASK_OPER_TUPLE                              \
  (OverlapOper, KineticOper, ElecPotOper, CartesianMultipoleOper<3u>, \
   CartesianMultipoleOper<3u>, CartesianMultipoleOper<3u>,            \
   SphericalMultipoleOper, σpVσpOper)
#define BOOST_PP_ONEBODY_TASK_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_ONEBODY_TASK_TUPLE)
#define BOOST_PP_ONEBODY_TASK_OPER_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_ONEBODY_TASK_OPER_TUPLE)

  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ONEBODY; ++d) {
#define BOOST_PP_ONEBODY_MCR1(r, data, elem) \
  taskmgr.add(task_label(BOOST_PP_STRINGIZE(elem), d));

    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR1, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR1
  }
#endif
#ifdef LIBINT_INCLUDE_ERI
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI; ++d) {
    taskmgr.add(task_label("eri", d));
  }
#endif

#ifdef LIBINT_INCLUDE_RKB_ERI
#define BOOST_PP_RKB_ERI_TASK_TUPLE \
  (coulomb_opop, opop_coulomb_opop, op_coulomb_op)
#define BOOST_PP_RKB_ERI_TASK_OPER_TUPLE \
  (CoulombσpσpOper, σpσpCoulombσpσpOper, opCoulombopOper)
#define BOOST_PP_RKB_ERI_TASK_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_RKB_ERI_TASK_TUPLE)
#define BOOST_PP_RKB_ERI_TASK_OPER_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_RKB_ERI_TASK_OPER_TUPLE)

  for (unsigned int d = 0; d <= LIBINT_INCLUDE_RKB_ERI; ++d) {
#define BOOST_PP_RKB_ERI_MCR1(r, data, elem) \
  taskmgr.add(task_label(BOOST_PP_STRINGIZE(elem), d));

    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR1, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR1
  }
#endif

#ifdef LIBINT_INCLUDE_ERI3
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI3; ++d) {
    taskmgr.add(task_label("3eri", d));
  }
#endif
#ifdef LIBINT_INCLUDE_ERI2
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI2; ++d) {
    taskmgr.add(task_label("2eri", d));
  }
#endif
#ifdef LIBINT_INCLUDE_G12
  taskmgr.add("r12kg12");
#if !LIBINT_USE_COMPOSITE_EVALUATORS
  taskmgr.add("r12_0_g12");
  taskmgr.add("r12_2_g12");
#endif
#endif
#ifdef LIBINT_INCLUDE_GENG12
  taskmgr.add("geng12");
#endif
#ifdef LIBINT_INCLUDE_G12DKH
  taskmgr.add("g12dkh");
#endif

  // use default parameters
  std::shared_ptr<CompilationParameters> cparams(new CompilationParameters);

  cparams->max_am("default", LIBINT_MAX_AM);
  cparams->max_am_opt("default", LIBINT_OPT_AM);
#if defined(LIBINT_MAX_AM_LIST)
  for (unsigned int d = 0; d < ntokens(LIBINT_MAX_AM_LIST, ','); ++d) {
    cparams->max_am(task_label("default", d),
                    token<unsigned int>(LIBINT_MAX_AM_LIST, ',', d));
  }
#endif
#if defined(LIBINT_OPT_AM_LIST)
  for (unsigned int d = 0; d < ntokens(LIBINT_OPT_AM_LIST, ','); ++d) {
    cparams->max_am_opt(task_label("default", d),
                        token<unsigned int>(LIBINT_OPT_AM_LIST, ',', d));
  }
#endif
  cparams->num_bf("default", 4);

#ifdef LIBINT_INCLUDE_ONEBODY
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ONEBODY; ++d) {
#if defined(LIBINT_ONEBODY_MAX_AM_LIST)
#define BOOST_PP_ONEBODY_MCR2(r, data, elem)   \
  cparams->max_am(                             \
      task_label(BOOST_PP_STRINGIZE(elem), d), \
                 token<unsigned int>(LIBINT_ONEBODY_MAX_AM_LIST, ',', d));
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR2, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR2
#elif defined(LIBINT_ONEBODY_MAX_AM)
#define BOOST_PP_ONEBODY_MCR3(r, data, elem)               \
  cparams->max_am(task_label(BOOST_PP_STRINGIZE(elem), d), \
                             LIBINT_ONEBODY_MAX_AM);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR3, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR3
#endif
#if defined(LIBINT_ONEBODY_OPT_AM_LIST)
#define BOOST_PP_ONEBODY_MCR4(r, data, elem)   \
  cparams->max_am_opt(                         \
      task_label(BOOST_PP_STRINGIZE(elem), d), \
                 token<unsigned int>(LIBINT_ONEBODY_OPT_AM_LIST, ',', d));
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR4, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR4
#elif defined(LIBINT_ONEBODY_OPT_AM)
#define BOOST_PP_ONEBODY_MCR5(r, data, elem)                   \
  cparams->max_am_opt(task_label(BOOST_PP_STRINGIZE(elem), d), \
                                 LIBINT_ONEBODY_OPT_AM);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR5, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR5
#endif
  }
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ONEBODY; ++d) {
#define BOOST_PP_ONEBODY_MCR6(r, data, elem) \
  cparams->num_bf(task_label(BOOST_PP_STRINGIZE(elem), d), 2);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_ONEBODY_MCR6, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR6
  }
#endif  // LIBINT_INCLUDE_ONEBODY

#ifdef LIBINT_INCLUDE_ERI
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI; ++d) {
#if defined(LIBINT_ERI_MAX_AM_LIST)
    cparams->max_am(task_label("eri", d),
                    token<unsigned int>(LIBINT_ERI_MAX_AM_LIST, ',', d));
#elif defined(LIBINT_ERI_MAX_AM)
    cparams->max_am(task_label("eri", d), LIBINT_ERI_MAX_AM);
#endif
#if defined(LIBINT_ERI_OPT_AM_LIST)
    cparams->max_am_opt(task_label("eri", d),
                        token<unsigned int>(LIBINT_ERI_OPT_AM_LIST, ',', d));
#elif defined(LIBINT_ERI_OPT_AM)
    cparams->max_am_opt(task_label("eri", d), LIBINT_ERI_OPT_AM);
#endif
  }
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI; ++d) {
    cparams->num_bf(task_label("eri", d), 4);
  }
#endif

#ifdef LIBINT_INCLUDE_RKB_ERI
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_RKB_ERI; ++d) {
#if defined(LIBINT_RKB_ERI_MAX_AM_LIST)
#define BOOST_PP_RKB_ERI_MCR2(r, data, elem)   \
  cparams->max_am(                             \
      task_label(BOOST_PP_STRINGIZE(elem), d), \
                 token<unsigned int>(LIBINT_RKB_ERI_MAX_AM_LIST, ',', d));
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR2, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR2
#elif defined(LIBINT_RKB_ERI_MAX_AM)
#define BOOST_PP_RKB_ERI_MCR3(r, data, elem)               \
  cparams->max_am(task_label(BOOST_PP_STRINGIZE(elem), d), \
                             LIBINT_RKB_ERI_MAX_AM);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR3, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR3
#endif
#if defined(LIBINT_RKB_ERI_OPT_AM_LIST)
#define BOOST_PP_RKB_ERI_MCR4(r, data, elem)   \
  cparams->max_am_opt(                         \
      task_label(BOOST_PP_STRINGIZE(elem), d), \
                 token<unsigned int>(LIBINT_RKB_ERI_OPT_AM_LIST, ',', d));
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR4, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR4
#elif defined(LIBINT_RKB_ERI_OPT_AM)
#define BOOST_PP_RKB_ERI_MCR5(r, data, elem)                   \
  cparams->max_am_opt(task_label(BOOST_PP_STRINGIZE(elem), d), \
                                 LIBINT_RKB_ERI_OPT_AM);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR5, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR5
#endif
  }
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_RKB_ERI; ++d) {
#define BOOST_PP_RKB_ERI_MCR6(r, data, elem) \
  cparams->num_bf(task_label(BOOST_PP_STRINGIZE(elem), d), 4);
    BOOST_PP_LIST_FOR_EACH(BOOST_PP_RKB_ERI_MCR6, _, BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR6
  }
#endif  // LIBINT_INCLUDE_RKB_ERI

#ifdef LIBINT_INCLUDE_ERI3
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI3; ++d) {
#if defined(LIBINT_ERI3_MAX_AM_LIST)
    cparams->max_am(task_label("3eri", d),
                    token<unsigned int>(LIBINT_ERI3_MAX_AM_LIST, ',', d));
#elif defined(LIBINT_ERI3_MAX_AM)
    cparams->max_am(task_label("3eri", d), LIBINT_ERI3_MAX_AM);
#endif
#if defined(LIBINT_ERI3_OPT_AM_LIST)
    cparams->max_am_opt(task_label("3eri", d),
                        token<unsigned int>(LIBINT_ERI3_OPT_AM_LIST, ',', d));
#elif defined(LIBINT_ERI3_OPT_AM)
    cparams->max_am_opt(task_label("3eri", d), LIBINT_ERI3_OPT_AM);
#endif

#if defined(LIBINT_MAX_AM_LIST)
    cparams->max_am(task_label("3eri", d),
                    cparams->max_am(task_label("default", d)), 1);
    cparams->max_am(task_label("3eri", d),
                    cparams->max_am(task_label("default", d)), 2);
#else
    cparams->max_am(task_label("3eri", d), cparams->max_am("default"), 1);
    cparams->max_am(task_label("3eri", d), cparams->max_am("default"), 2);
#endif
  }
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI3; ++d) {
    cparams->num_bf(task_label("3eri", d), 3);
  }
#endif
#ifdef LIBINT_INCLUDE_ERI2
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI2; ++d) {
#if defined(LIBINT_ERI2_MAX_AM_LIST)
    cparams->max_am(task_label("2eri", d),
                    token<unsigned int>(LIBINT_ERI2_MAX_AM_LIST, ',', d));
#elif defined(LIBINT_ERI2_MAX_AM)
    cparams->max_am(task_label("2eri", d), LIBINT_ERI2_MAX_AM);
#endif
#if defined(LIBINT_ERI2_OPT_AM_LIST)
    cparams->max_am_opt(task_label("2eri", d),
                        token<unsigned int>(LIBINT_ERI2_OPT_AM_LIST, ',', d));
#elif defined(LIBINT_ERI2_OPT_AM)
    cparams->max_am_opt(task_label("2eri", d), LIBINT_ERI2_OPT_AM);
#endif
  }
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI2; ++d) {
    cparams->num_bf(task_label("2eri", d), 2);
  }
#endif
#ifdef LIBINT_INCLUDE_G12
#ifndef LIBINT_G12_MAX_AM
#define LIBINT_MAX_AM LIBINT_G12_MAX_AM
#endif
#ifndef LIBINT_G12_OPT_AM
#define LIBINT_OPT_AM LIBINT_G12_OPT_AM
#endif
  cparams->max_am("r12kg12", LIBINT_G12_MAX_AM);
  cparams->max_am_opt("r12kg12", LIBINT_G12_OPT_AM);
  cparams->num_bf("r12kg12", 4);
#if !LIBINT_USE_COMPOSITE_EVALUATORS
  cparams->max_am("r12_0_g12", LIBINT_G12_MAX_AM);
  cparams->max_am_opt("r12_0_g12", LIBINT_G12_OPT_AM);
  cparams->num_bf("r12_0_g12", 4);
  cparams->max_am("r12_2_g12", LIBINT_G12_MAX_AM);
  cparams->max_am_opt("r12_2_g12", LIBINT_G12_OPT_AM);
  cparams->num_bf("r12_2_g12", 4);
#endif
#endif
#ifdef LIBINT_INCLUDE_G12DKH
#ifndef LIBINT_G12DKH_MAX_AM
#define LIBINT_MAX_AM LIBINT_G12DKH_MAX_AM
#endif
#ifndef LIBINT_G12DKH_OPT_AM
#define LIBINT_OPT_AM LIBINT_G12DKH_OPT_AM
#endif
  cparams->max_am("g12dkh", LIBINT_G12DKH_MAX_AM);
  cparams->max_am_opt("g12dkh", LIBINT_G12DKH_OPT_AM);
#endif
#if LIBINT_ENABLE_UNROLLING
  cparams->unroll_threshold(LIBINT_ENABLE_UNROLLING);
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
#ifdef LIBINT_ALIGN_SIZE
  cparams->align_size(LIBINT_ALIGN_SIZE);
#endif
#if LIBINT_FLOP_COUNT
  cparams->count_flops(true);
#endif
#if LIBINT_PROFILE
  cparams->profile(true);
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
#ifdef LIBINT_USER_DEFINED_REAL
  {
    const std::string realtype(LIBINT_USER_DEFINED_REAL);
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
  { cparams->single_evaltype(true); }
#else
  throw std::runtime_error("Cannot generate specialized evaluator types yet");
#endif

  // initialize code context to produce library API
  std::shared_ptr<CodeContext> icontext(new CppCodeContext(cparams));
  // initialize object to generate interface
  std::shared_ptr<Libint2Iface> iface(new Libint2Iface(cparams, icontext));

  print_header(os);
  print_config(os);
  // transfer some configuration parameters to the generated library API
  iface->to_params(iface->macro_define("USE_COMPOSITE_EVALUATORS",
                                       LIBINT_USE_COMPOSITE_EVALUATORS));
  iface->to_params(
      iface->macro_define("CARTGAUSS_MAX_AM", LIBINT_CARTGAUSS_MAX_AM));
  iface->to_params(
      iface->macro_define("CGSHELL_ORDERING", LIBINT_CGSHELL_ORDERING));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_STANDARD",
                                       LIBINT_CGSHELL_ORDERING_STANDARD));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_INTV3",
                                       LIBINT_CGSHELL_ORDERING_INTV3));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_GAMESS",
                                       LIBINT_CGSHELL_ORDERING_GAMESS));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_ORCA",
                                       LIBINT_CGSHELL_ORDERING_ORCA));
  iface->to_params(iface->macro_define("CGSHELL_ORDERING_BAGEL",
                                       LIBINT_CGSHELL_ORDERING_BAGEL));
  iface->to_params(iface->macro_define("SHELLQUARTET_SET", LIBINT_SHELL_SET));
  iface->to_params(iface->macro_define("SHELLQUARTET_SET_STANDARD",
                                       LIBINT_SHELL_SET_STANDARD));
  iface->to_params(
      iface->macro_define("SHELLQUARTET_SET_ORCA", LIBINT_SHELL_SET_ORCA));
#if defined(LIBINT_MAX_AM_LIST)
  for (unsigned int d = 0; d < ntokens(LIBINT_MAX_AM_LIST, ','); ++d) {
    {
      std::ostringstream oss;
      oss << "MAX_AM";
      if (d > 0) oss << d;
      iface->to_params(iface->macro_define(
          oss.str(), token<unsigned int>(LIBINT_MAX_AM_LIST, ',', d)));
    }
    {
      std::ostringstream oss;
      oss << "MAX_AM_default";
      if (d > 0) oss << d;
      iface->to_params(iface->macro_define(
          oss.str(), token<unsigned int>(LIBINT_MAX_AM_LIST, ',', d)));
    }
  }
#else
  iface->to_params(iface->macro_define("MAX_AM", LIBINT_MAX_AM));
  int max_deriv = 0;
#ifdef LIBINT_INCLUDE_ONEBODY
  max_deriv = std::max(LIBINT_INCLUDE_ONEBODY, max_deriv);
#endif
#ifdef LIBINT_INCLUDE_ERI
  max_deriv = std::max(LIBINT_INCLUDE_ERI, max_deriv);
#endif
#ifdef LIBINT_INCLUDE_RKB_ERI
  max_deriv = std::max(LIBINT_INCLUDE_RKB_ERI, max_deriv);
#endif
#ifdef LIBINT_INCLUDE_ERI3
  max_deriv = std::max(LIBINT_INCLUDE_ERI3, max_deriv);
#endif
#ifdef LIBINT_INCLUDE_ERI2
  max_deriv = std::max(LIBINT_INCLUDE_ERI2, max_deriv);
#endif
  for (int d = 0; d <= max_deriv; ++d) {
    std::ostringstream oss;
    oss << "MAX_AM_default";
    if (d > 0) oss << d;
    iface->to_params(iface->macro_define(oss.str(), LIBINT_MAX_AM));
  }
#endif
  cparams->print(os);

  g_progress.start();

#ifdef LIBINT_INCLUDE_ONEBODY
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ONEBODY; ++d) {
#define BOOST_PP_ONEBODY_MCR7(r, data, i, elem)                              \
  build_onebody_1b_1k<BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_TASK_OPER_LIST, i)>( \
      os, BOOST_PP_STRINGIZE(elem), cparams, iface, d);
    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_ONEBODY_MCR7, _,
                             BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_ONEBODY_MCR7
  }
#endif
#ifdef LIBINT_INCLUDE_ERI
#if !USE_GENERIC_ERI_BUILD
  build_TwoPRep_2b_2k<TwoPRep>(os, "eri", cparams, iface);
#else
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI; ++d) {
    build_TwoPRep_2b_2k<TwoPRep>(os, "eri", cparams, iface, d);
  }
#endif
#endif

#ifdef LIBINT_INCLUDE_RKB_ERI
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_RKB_ERI; ++d) {
#define BOOST_PP_RKB_ERI_MCR7(r, data, i, elem)                              \
  build_TwoPRep_2b_2k<BOOST_PP_LIST_AT(BOOST_PP_RKB_ERI_TASK_OPER_LIST, i)>( \
      os, BOOST_PP_STRINGIZE(elem), cparams, iface, d);
    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_RKB_ERI_MCR7, _,
                             BOOST_PP_RKB_ERI_TASK_LIST)
#undef BOOST_PP_RKB_ERI_MCR7
  }
#endif

#ifdef LIBINT_INCLUDE_ERI3
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI3; ++d) {
    build_TwoPRep_1b_2k(os, cparams, iface, d);
  }
#if LIBINT_ERI3_PURE_SH
  iface->to_params(iface->macro_define("ERI3_PURE_SH", 1));
#endif
#endif
#ifdef LIBINT_INCLUDE_ERI2
  for (unsigned int d = 0; d <= LIBINT_INCLUDE_ERI2; ++d) {
    build_TwoPRep_1b_1k(os, cparams, iface, d);
  }
#if LIBINT_ERI2_PURE_SH
  iface->to_params(iface->macro_define("ERI2_PURE_SH", 1));
#endif
#endif
#ifdef LIBINT_INCLUDE_G12
#if LIBINT_USE_COMPOSITE_EVALUATORS
  build_R12kG12_2b_2k(os, cparams, iface);
#else
  build_R12kG12_2b_2k_separate(os, cparams, iface);
#endif
#endif
#ifdef LIBINT_INCLUDE_G12DKH
  build_G12DKH_2b_2k(os, cparams, iface);
#endif

  g_progress.finish();

  // Generate code for the set-level RRs
  std::deque<std::string> decl_filenames, def_filenames;
  generate_rr_code(os, cparams, decl_filenames, def_filenames);

#if DEBUG
  // print out the external symbols found for each task
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  for (tciter t = taskmgr.first(); t != tend; ++t) {
    const std::shared_ptr<TaskExternSymbols> tsymbols = t->symbols();
    typedef TaskExternSymbols::SymbolList SymbolList;
    const SymbolList& symbols = tsymbols->symbols();
    // print out the labels
    std::cout << "Recovered labels for task " << t->label() << std::endl;
    typedef SymbolList::const_iterator citer;
    citer end = symbols.end();
    for (citer s = symbols.begin(); s != end; ++s) std::cout << *s << std::endl;
  }
#endif

  // transfer some library configuration to library API
  config_to_api(cparams, iface);

  os << "Compilation finished. Goodbye." << endl;
}

void print_header(std::ostream& os) {
  os << "----------------------------------------------" << endl;
  os << " build_libint2: optimizing integrals compiler " << endl;
  os << "                        by Edward F. Valeev   " << endl;
  os << "                and ideas by Justin Fermann   " << endl;
  os << "                         and Curtis Janssen   " << endl;
  os << "----------------------------------------------" << endl << endl;
}

void print_config(std::ostream& os) {
#ifdef LIBINT_INCLUDE_ONEBODY
  os << "Will support one-body integrals";
  if (LIBINT_INCLUDE_ONEBODY > 0)
    os << "(deriv order = " << LIBINT_INCLUDE_ONEBODY << ")";
  os << endl;
#endif
#ifdef LIBINT_INCLUDE_ERI
  os << "Will support ERI";
  if (LIBINT_INCLUDE_ERI > 0)
    os << "(deriv order = " << LIBINT_INCLUDE_ERI << ")";
  os << endl;
#endif
#ifdef LIBINT_INCLUDE_ERI3
  os << "Will support 3-center ERI";
  if (LIBINT_INCLUDE_ERI3 > 0)
    os << "(deriv order = " << LIBINT_INCLUDE_ERI3 << ")";
  os << endl;
#endif
#ifdef LIBINT_INCLUDE_ERI2
  os << "Will support 2-center ERI";
  if (LIBINT_INCLUDE_ERI2 > 0)
    os << "(deriv order = " << LIBINT_INCLUDE_ERI2 << ")";
  os << endl;
#endif
#ifdef LIBINT_INCLUDE_G12
  os << "Will support G12 (commutators = ";
#if LIBINT_SUPPORT_T1G12
  os << "yes";
#else
  os << "false";
#endif
  os << ")" << endl;
#endif
#ifdef LIBINT_INCLUDE_G12DKH
  os << "Will support G12DKH" << endl;
#endif
#ifdef LIBINT_INCLUDE_RKB_ERI
  os << "Will support restricted kinetically balance (RKB) 4-center ERIs "
     << std::endl;
  if (LIBINT_INCLUDE_RKB_ERI > 0)
    os << "(deriv order = " << LIBINT_INCLUDE_RKB_ERI << ")";
  os << endl;
#endif
}

#if defined(LIBINT_INCLUDE_ERI) || defined(LIBINT_INCLUDE_RKB_ERI)
template <typename OperType>
static void build_TwoPRep_2b_2k(
    std::ostream& os, std::string label,
    const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface, unsigned int deriv_level) {
  typedef GenIntegralSet_11_11<CGShell, OperType, mType> TwoBody_sh_11_11;
  typedef typename OperType::Descriptor OperDescrType;

  const std::string task = task_label(label, deriv_level);

  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define(std::string("MAX_AM_") + task, lmax));

  const auto nullaux = typename TwoBody_sh_11_11::AuxIndexType(0u);
  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  std::shared_ptr<DirectedGraph> dg_xxxx(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy());
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

  // opCoulombop has only a 2-fold bra↔ket-swap symmetry (with (a,b)↔(b,a)
  // component remap). Within-side particle swap is NOT a symmetry because σ·p
  // attaches to one specific function per side (ν on bra, λ on ket); swapping
  // moves the operator to a different physical center that IBP cannot recover
  // when centers differ. Emit code for every (la,lb,lc,ld) combination to
  // avoid triggering any within-side swap at runtime.
  bool p1_p2_swappable = !std::is_same<OperType, CoulombσpσpOper>::value &&
                         !std::is_same<OperType, opCoulombopOper>::value;
  bool bra_ket_coswappable = std::is_same<OperType, σpσpCoulombσpσpOper>::value;

  // Note: la, lb, lc, ld generate code for chemist notation (ab|O|cd), where O
  // is a two-body operator.
  for (unsigned int la = 0; la <= lmax; la++) {
    for (unsigned int lb = 0; lb <= lmax; lb++) {
      for (unsigned int lc = 0; lc <= lmax; lc++) {
        for (unsigned int ld = 0; ld <= lmax; ld++) {
          // opCoulombop has only a bra↔ket (particle 1↔2) swap symmetry;
          // within-side swap is NOT a symmetry (σ·p would move to the wrong
          // physical center). Canonical form: la+lb <= lc+ld only
          // (ORCA: la+lb >= lc+ld). Use a dedicated predicate so within-side
          // orderings are not reduced away.
          if constexpr (std::is_same<OperType, opCoulombopOper>::value) {
#if LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD
            if (!(la + lb <= lc + ld)) continue;
#else
            if (!(la + lb >= lc + ld)) continue;
#endif
          } else {
            if (!ShellQuartetSetPredicate<static_cast<ShellSetType>(
                    LIBINT_SHELL_SET)>::value(la, lb, lc, ld, p1_p2_swappable,
                                              bra_ket_coswappable))
              continue;
          }

          // std::shared_ptr<Tactic> tactic(new ParticleDirectionTactic(la+lb >
          // lc+ld ? false : true));
          std::shared_ptr<Tactic> tactic(
              new FourCenter_OS_Tactic(la, lb, lc, ld));

#if STUDY_MEMORY_USAGE
          const int lim = 1;
          if (!(la == lim && lb == lim && lc == lim && ld == lim)) continue;
#endif
          // this will hold all target shell sets
          std::vector<std::shared_ptr<TwoBody_sh_11_11>> targets;

          /////////////////////////////////
          // loop over operator components
          /////////////////////////////////
          std::vector<OperDescrType> descrs(1);
          if constexpr (std::is_same<OperType, CoulombσpσpOper>::value) {
            // reset descriptors array
            descrs.resize(0);
            // iterate over 4 quaternion components (single spin space)
            for (int p = 0; p != 4; ++p) {
              descrs.emplace_back(OperDescrType(p));
            }
          }
          if constexpr (std::is_same<OperType, σpσpCoulombσpσpOper>::value) {
            // reset descriptors array
            descrs.resize(0);
            // iterate over 16 components (tensor product of two spin spaces)
            for (int p = 0; p != 16; ++p) {
              descrs.emplace_back(OperDescrType(p));
            }
          }
          if constexpr (std::is_same<OperType, opCoulombopOper>::value) {
            // reset descriptors array
            descrs.resize(0);
            // iterate over 9 components (3x3 Cartesian dyadic: bra-dir ×
            // ket-dir)
            for (int p = 0; p != 9; ++p) {
              descrs.emplace_back(OperDescrType(p));
            }
          }

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la, lb), max(lc, ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const auto nopers = descrs.size();
          const bool need_to_unroll =
              nopers * l_to_cgshellsize(la) * l_to_cgshellsize(lb) *
                  l_to_cgshellsize(lc) * l_to_cgshellsize(ld) <=
              cparams->unroll_threshold();
          const unsigned int unroll_threshold =
              need_to_optimize && need_to_unroll
                  ? std::numeric_limits<unsigned int>::max()
                  : 0;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          // For multi-component operators (RKB), components share no
          // intermediates, so CSE/condense_expr is pure overhead — disable.
          const bool do_optimize = (nopers > 1) ? false : need_to_optimize;
          dg_xxxx->registry()->do_cse(do_optimize);
          dg_xxxx->registry()->condense_expr(
              do_optimize ? condense_expr(cparams->unroll_threshold(),
                                          cparams->max_vector_length() > 1)
                          : false);
          //  Need to accumulate integrals?
          dg_xxxx->registry()->accumulate_targets(
              cparams->accumulate_targets());
          // need to profile?
          if (cparams->profile()) {
            dg_xxxx->registry()->current_timer(0);
          }

          ////////////
          // loop over unique derivative index combinations
          ////////////
          // NB translational invariance is now handled by CR_DerivGauss
          CartesianDerivIterator<4> diter(deriv_level);

          bool last_deriv = false;
          do {
            CGShell a(la);
            CGShell b(lb);
            CGShell c(lc);
            CGShell d(ld);

            for (unsigned int i = 0; i < 4; ++i) {
              for (unsigned int xyz = 0; xyz < 3; ++xyz) {
                if (i == 0) a.deriv().inc(xyz, (*diter).at(3 * i + xyz));
                if (i == 1) b.deriv().inc(xyz, (*diter).at(3 * i + xyz));
                if (i == 2) c.deriv().inc(xyz, (*diter).at(3 * i + xyz));
                if (i == 3) d.deriv().inc(xyz, (*diter).at(3 * i + xyz));
              }
            }

            // operator component loop
            for (unsigned int op = 0; op != descrs.size(); ++op) {
              OperType oper(descrs[op]);

              std::shared_ptr<TwoBody_sh_11_11> abcd =
                  TwoBody_sh_11_11::Instance(a, b, c, d, nullaux, oper);
              targets.push_back(abcd);
            }

            last_deriv = diter.last();
            if (!last_deriv) diter.next();
          } while (!last_deriv);
          // make label that characterizes this set of targets
          // use the label of the nondifferentiated integral as a base
          std::string abcd_label;
          {
            CGShell a(la);
            CGShell b(lb);
            CGShell c(lc);
            CGShell d(ld);

            if constexpr (std::is_same<OperType, TwoPRep>::value) {
              OperType oper;
              oper = OperType(descrs[0]);
              std::shared_ptr<TwoBody_sh_11_11> abcd =
                  TwoBody_sh_11_11::Instance(a, b, c, d, nullaux, oper);
              abcd_label = abcd->label();
            } else {
              std::ostringstream oss;
              oss << "_" << a.label() << "_" << b.label();
              oss << "_" << label;
              oss << "_" << c.label() << "_" << d.label();
              abcd_label = oss.str();
            }
          }
          // + derivative level (if deriv_level > 0)
          std::string eval_label;
          {
            eval_label = "";
            if (deriv_level != 0) {
              std::ostringstream oss;
              oss << "deriv" << deriv_level;
              eval_label += oss.str();
            }
            eval_label += abcd_label;
          }

          g_progress.current_task = eval_label;
          g_progress.print();

          std::string prefix(cparams->source_directory());
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // append all targets to the graph
          for (auto it = targets.begin(); it != targets.end(); ++it) {
            std::shared_ptr<DGVertex> t_ptr =
                std::dynamic_pointer_cast<DGVertex, TwoBody_sh_11_11>(*it);
            dg_xxxx->append_target(t_ptr);
          }

          // this will generate code for these targets, and potentially
          // generate code for its prerequisites
          GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames, prefix, eval_label,
                       false);

          // update max stack size and # of targets
          const std::shared_ptr<TaskParameters>& tparams =
              taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(targets.size());
          // os << " Max memory used = " << memman->max_memory_used() <<
          // std::endl;

          // set pointer to the top-level evaluator function
          ostringstream oss;
          oss << context->label_to_function_name("libint2_build_" + task) << "["
              << la << "][" << lb << "][" << lc << "][" << ld
              << "] = " << context->label_to_function_name(eval_label)
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for (auto& decl_filename : decl_filenames) {
            oss.str("");
            oss << "#include <" << decl_filename << ">" << endl;
            iface->to_int_iface(oss.str());
          }

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxxx->reset();
          memman->reset();

          ++g_progress.done;
          g_progress.print();

        }  // end of d loop
      }    // end of c loop
    }      // end of b loop
  }        // end of a loop
}

#endif  // LIBINT_INCLUDE_ERI || LIBINT_INCLUDE_RKB_ERI

#ifdef LIBINT_INCLUDE_ERI3

void build_TwoPRep_1b_2k(std::ostream& os,
                         const std::shared_ptr<CompilationParameters>& cparams,
                         std::shared_ptr<Libint2Iface>& iface,
                         unsigned int deriv_level) {
  const std::string task = task_label("3eri", deriv_level);
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  const unsigned int lmax = cparams->max_am(task);
  const unsigned int lmax_default =
      const_cast<const CompilationParameters*>(cparams.get())->max_am(task, 1);
  if (lmax != lmax_default)
    iface->to_params(
        iface->macro_define(std::string("CENTER_DEPENDENT_MAX_AM_") + task, 1));

  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define(std::string("MAX_AM_") + task, lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  std::shared_ptr<DirectedGraph> dg_xxx(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy());
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

  for (unsigned int lbra = 0; lbra <= lmax; lbra++) {
    for (unsigned int lc = 0; lc <= lmax_default; lc++) {
      for (unsigned int ld = 0; ld <= lmax_default; ld++) {
        // eliminate some cases depending on the desired convention
        if (!ShellTripletSetPredicate<static_cast<ShellSetType>(
                LIBINT_SHELL_SET)>::value(lbra, lc, ld))
          continue;

        // I will use 4-center recurrence relations and integrals, and have one
        // center carry an s function unfortunately, depending on the direction
        // in which the build goes it must be A(0) or B(1)
        const unsigned int dummy_center =
            (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 0 : 1;

        // std::shared_ptr<Tactic> tactic(new ParticleDirectionTactic(lbra >
        // lc+ld ? false : true));
        std::shared_ptr<Tactic> tactic(
            new FourCenter_OS_Tactic(dummy_center == 0 ? 0 : lbra,
                                     dummy_center == 1 ? 0 : lbra, lc, ld));

#if STUDY_MEMORY_USAGE
        const int lim = 1;
        if (!(lbra == lim && lc == lim && ld == lim)) continue;
#endif

        // unroll only if max_am <= cparams->max_am_opt(task)
        using std::max;
        const unsigned int max_am = max(max(lc, ld), lbra);
        const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
        const bool need_to_unroll = l_to_cgshellsize(lbra) *
                                        l_to_cgshellsize(lc) *
                                        l_to_cgshellsize(ld) <=
                                    cparams->unroll_threshold();
        const unsigned int unroll_threshold =
            need_to_optimize && need_to_unroll
                ? std::numeric_limits<unsigned int>::max()
                : 0;
        dg_xxx->registry()->unroll_threshold(unroll_threshold);
        dg_xxx->registry()->do_cse(need_to_optimize);
        dg_xxx->registry()->condense_expr(condense_expr(
            cparams->unroll_threshold(), cparams->max_vector_length() > 1));
        // dg_xxx->registry()->condense_expr(true);
        //  Need to accumulate integrals?
        dg_xxx->registry()->accumulate_targets(cparams->accumulate_targets());

        ////////////
        // loop over unique derivative index combinations
        ////////////
        // NB translational invariance is now handled by CR_DerivGauss
        CartesianDerivIterator<3> diter(deriv_level);
        std::vector<std::shared_ptr<TwoPRep_sh_11_11>> targets;
        bool last_deriv = false;
        do {
          CGShell a = (dummy_center == 0) ? CGShell::unit() : CGShell(lbra);
          CGShell b = (dummy_center == 1) ? CGShell::unit() : CGShell(lbra);
          CGShell c(lc);
          CGShell d(ld);
#if LIBINT_ERI3_PURE_SH
          if (dummy_center == 1 && deriv_level == 0) a.pure_sh(true);
          if (dummy_center == 0 && deriv_level == 0) b.pure_sh(true);
#endif

          unsigned int center = 0;
          for (unsigned int i = 0; i < 4; ++i) {
            if (i == dummy_center) continue;
            for (unsigned int xyz = 0; xyz < 3; ++xyz) {
              if (i == 0) a.deriv().inc(xyz, (*diter).at(3 * center + xyz));
              if (i == 1) b.deriv().inc(xyz, (*diter).at(3 * center + xyz));
              if (i == 2) c.deriv().inc(xyz, (*diter).at(3 * center + xyz));
              if (i == 3) d.deriv().inc(xyz, (*diter).at(3 * center + xyz));
            }
            ++center;
          }

          // use 4-center integrals
          std::shared_ptr<TwoPRep_sh_11_11> abcd =
              TwoPRep_sh_11_11::Instance(a, b, c, d, mType(0u));
          targets.push_back(abcd);
          last_deriv = diter.last();
          if (!last_deriv) diter.next();
        } while (!last_deriv);
        // append all derivatives as targets to the graph
        for (std::vector<std::shared_ptr<TwoPRep_sh_11_11>>::const_iterator t =
                 targets.begin();
             t != targets.end(); ++t) {
          std::shared_ptr<DGVertex> t_ptr =
              std::dynamic_pointer_cast<DGVertex, TwoPRep_sh_11_11>(*t);
          dg_xxx->append_target(t_ptr);
        }

        // make label that characterizes this set of targets
        // use the label of the nondifferentiated integral as a base
        std::string abcd_label;
        {
          CGShell a = (dummy_center == 0) ? CGShell::unit() : CGShell(lbra);
          CGShell b = (dummy_center == 1) ? CGShell::unit() : CGShell(lbra);
          CGShell c(lc);
          CGShell d(ld);
#if LIBINT_ERI3_PURE_SH
          if (dummy_center == 1 && deriv_level == 0) a.pure_sh(true);
          if (dummy_center == 0 && deriv_level == 0) b.pure_sh(true);
#endif
          std::shared_ptr<TwoPRep_sh_11_11> abcd =
              TwoPRep_sh_11_11::Instance(a, b, c, d, mType(0u));
          abcd_label = abcd->label();
        }
        // + derivative level (if deriv_level > 0)
        std::string label;
        {
          label = "";
          if (deriv_level != 0) {
            std::ostringstream oss;
            oss << "deriv" << deriv_level;
            label += oss.str();
          }
          label += "eri3";
          label += abcd_label;
        }

        g_progress.current_task = label;
        g_progress.print();

        std::string prefix(cparams->source_directory());
        std::deque<std::string> decl_filenames;
        std::deque<std::string> def_filenames;

        // this will generate code for this targets, and potentially generate
        // code for its prerequisites
        GenerateCode(dg_xxx, context, cparams, strat, tactic, memman,
                     decl_filenames, def_filenames, prefix, label, false);

        // update max stack size and # of targets
        const std::shared_ptr<TaskParameters>& tparams =
            taskmgr.current().params();
        tparams->max_stack_size(max_am, memman->max_memory_used());
        tparams->max_ntarget(targets.size());
        // os << " Max memory used = " << memman->max_memory_used() <<
        // std::endl;

        // set pointer to the top-level evaluator function
        ostringstream oss;
        oss << context->label_to_function_name("libint2_build_" + task) << "["
            << lbra << "][" << lc << "][" << ld
            << "] = " << context->label_to_function_name(label)
            << context->end_of_stat() << endl;
        iface->to_static_init(oss.str());

        // need to declare this function internally
        for (auto& decl_filename : decl_filenames) {
          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());
        }

#if DEBUG
        os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
        dg_xxx->reset();
        memman->reset();
        ++g_progress.done;
        g_progress.print();
      }  // end of d loop
    }    // end of c loop
  }      // end of bra loop
}
#endif  // LIBINT_INCLUDE_ERI3

#ifdef LIBINT_INCLUDE_ERI2

void build_TwoPRep_1b_1k(std::ostream& os,
                         const std::shared_ptr<CompilationParameters>& cparams,
                         std::shared_ptr<Libint2Iface>& iface,
                         unsigned int deriv_level) {
  const std::string task = task_label("2eri", deriv_level);
  typedef TwoPRep_11_11_sq TwoPRep_sh_11_11;
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define(std::string("MAX_AM_") + task, lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  std::shared_ptr<DirectedGraph> dg_xxx(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy());
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

  for (unsigned int lbra = 0; lbra <= lmax; lbra++) {
    for (unsigned int lket = 0; lket <= lmax; lket++) {
      // I will use 4-center recurrence relations and integrals, and have two
      // centers carry an s function unfortunately, depending on the direction
      // in which the build goes it must be A(0) and C(2) or B(1) and D(3)
      const unsigned int dummy_center1 =
          (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 0 : 1;
      const unsigned int dummy_center2 =
          (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? 2 : 3;

      // std::shared_ptr<Tactic> tactic(new ParticleDirectionTactic(lbra > lket
      // ? false : true));
      std::shared_ptr<Tactic> tactic(new FourCenter_OS_Tactic(
          dummy_center1 == 0 ? 0 : lbra, dummy_center1 == 1 ? 0 : lbra,
          dummy_center2 == 2 ? 0 : lket, dummy_center2 == 3 ? 0 : lket));

#if STUDY_MEMORY_USAGE
      const int lim = 1;
      if (!(lbra == lim && lket == lim)) continue;
#endif

      // unroll only if max_am <= cparams->max_am_opt(task)
      using std::max;
      const unsigned int max_am = max(lbra, lket);
      const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
      const bool need_to_unroll =
          l_to_cgshellsize(lbra) * l_to_cgshellsize(lket) <=
          cparams->unroll_threshold();
      const unsigned int unroll_threshold =
          need_to_optimize && need_to_unroll
              ? std::numeric_limits<unsigned int>::max()
              : 0;
      dg_xxx->registry()->unroll_threshold(unroll_threshold);
      dg_xxx->registry()->do_cse(need_to_optimize);
      dg_xxx->registry()->condense_expr(condense_expr(
          cparams->unroll_threshold(), cparams->max_vector_length() > 1));
      // Need to accumulate integrals?
      dg_xxx->registry()->accumulate_targets(cparams->accumulate_targets());

      ////////////
      // loop over unique derivative index combinations
      ////////////
      // NB translational invariance is now handled by CR_DerivGauss
      CartesianDerivIterator<2> diter(deriv_level);
      std::vector<std::shared_ptr<TwoPRep_sh_11_11>> targets;
      bool last_deriv = false;
      do {
        CGShell a = (dummy_center1 == 0) ? CGShell::unit() : CGShell(lbra);
        CGShell b = (dummy_center1 == 1) ? CGShell::unit() : CGShell(lbra);
        CGShell c = (dummy_center2 == 2) ? CGShell::unit() : CGShell(lket);
        CGShell d = (dummy_center2 == 3) ? CGShell::unit() : CGShell(lket);
#if LIBINT_ERI2_PURE_SH
        if (dummy_center1 == 1 && deriv_level == 0) a.pure_sh(true);
        if (dummy_center1 == 0 && deriv_level == 0) b.pure_sh(true);
        if (dummy_center2 == 3 && deriv_level == 0) c.pure_sh(true);
        if (dummy_center2 == 2 && deriv_level == 0) d.pure_sh(true);
#endif

        unsigned int center = 0;
        for (unsigned int i = 0; i < 4; ++i) {
          if (i == dummy_center1 || i == dummy_center2) continue;
          for (unsigned int xyz = 0; xyz < 3; ++xyz) {
            if (i == 0) a.deriv().inc(xyz, (*diter).at(3 * center + xyz));
            if (i == 1) b.deriv().inc(xyz, (*diter).at(3 * center + xyz));
            if (i == 2) c.deriv().inc(xyz, (*diter).at(3 * center + xyz));
            if (i == 3) d.deriv().inc(xyz, (*diter).at(3 * center + xyz));
          }
          ++center;
        }

        // use 4-center integrals
        std::shared_ptr<TwoPRep_sh_11_11> abcd =
            TwoPRep_sh_11_11::Instance(a, b, c, d, mType(0u));
        targets.push_back(abcd);
        last_deriv = diter.last();
        if (!last_deriv) diter.next();
      } while (!last_deriv);
      // append all derivatives as targets to the graph
      for (std::vector<std::shared_ptr<TwoPRep_sh_11_11>>::const_iterator t =
               targets.begin();
           t != targets.end(); ++t) {
        std::shared_ptr<DGVertex> t_ptr =
            std::dynamic_pointer_cast<DGVertex, TwoPRep_sh_11_11>(*t);
        dg_xxx->append_target(t_ptr);
      }

      // make label that characterizes this set of targets
      // use the label of the nondifferentiated integral as a base
      std::string abcd_label;
      {
        CGShell a = (dummy_center1 == 0) ? CGShell::unit() : CGShell(lbra);
        CGShell b = (dummy_center1 == 1) ? CGShell::unit() : CGShell(lbra);
        CGShell c = (dummy_center2 == 2) ? CGShell::unit() : CGShell(lket);
        CGShell d = (dummy_center2 == 3) ? CGShell::unit() : CGShell(lket);
#if LIBINT_ERI2_PURE_SH
        if (dummy_center1 == 1 && deriv_level == 0) a.pure_sh(true);
        if (dummy_center1 == 0 && deriv_level == 0) b.pure_sh(true);
        if (dummy_center2 == 3 && deriv_level == 0) c.pure_sh(true);
        if (dummy_center2 == 2 && deriv_level == 0) d.pure_sh(true);
#endif
        std::shared_ptr<TwoPRep_sh_11_11> abcd =
            TwoPRep_sh_11_11::Instance(a, b, c, d, mType(0u));
        abcd_label = abcd->label();
      }
      // + derivative level (if deriv_level > 0)
      std::string label;
      {
        label = "";
        if (deriv_level != 0) {
          std::ostringstream oss;
          oss << "deriv" << deriv_level;
          label += oss.str();
        }
        label += "eri2";
        label += abcd_label;
      }

      g_progress.current_task = label;
      g_progress.print();
      std::cout.flush();

      std::string prefix(cparams->source_directory());
      std::deque<std::string> decl_filenames;
      std::deque<std::string> def_filenames;

      // this will generate code for this targets, and potentially generate
      // code for its prerequisites
      GenerateCode(dg_xxx, context, cparams, strat, tactic, memman,
                   decl_filenames, def_filenames, prefix, label, false);

      // update max stack size and # of targets
      const std::shared_ptr<TaskParameters>& tparams =
          taskmgr.current().params();
      tparams->max_stack_size(max_am, memman->max_memory_used());
      tparams->max_ntarget(targets.size());
      // os << " Max memory used = " << memman->max_memory_used() <<
      // std::endl;

      // set pointer to the top-level evaluator function
      ostringstream oss;
      oss << context->label_to_function_name("libint2_build_" + task) << "["
          << lbra << "][" << lket
          << "] = " << context->label_to_function_name(label)
          << context->end_of_stat() << endl;
      iface->to_static_init(oss.str());

      // need to declare this function internally
      for (auto& decl_filename : decl_filenames) {
        oss.str("");
        oss << "#include <" << decl_filename << ">" << endl;
        iface->to_int_iface(oss.str());
      }

#if DEBUG
      os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
      dg_xxx->reset();
      memman->reset();
    }  // end of ket loop
  }    // end of bra loop
}
#endif  // LIBINT_INCLUDE_ERI2

#ifdef LIBINT_INCLUDE_G12
void build_R12kG12_2b_2k(std::ostream& os,
                         const std::shared_ptr<CompilationParameters>& cparams,
                         std::shared_ptr<Libint2Iface>& iface) {
  const std::string task("r12kg12");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_R12kG12", lmax));
#if LIBINT_SUPPORT_T1G12
  iface->to_params(iface->macro_define("SUPPORT_T1G12", 1));
#else
  iface->to_params(iface->macro_define("SUPPORT_T1G12", 0));
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
  std::shared_ptr<DirectedGraph> dg_xxxx(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy);
  std::shared_ptr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  // std::shared_ptr<Tactic> tactic(new RandomChoiceTactic());
  // std::shared_ptr<Tactic> tactic(new FewestNewVerticesTactic(dg_xxxx));
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

  for (unsigned int la = 0; la <= lmax; la++) {
    for (unsigned int lb = 0; lb <= lmax; lb++) {
      for (unsigned int lc = 0; lc <= lmax; lc++) {
        for (unsigned int ld = 0; ld <= lmax; ld++) {
          if (la < lb || lc < ld || la + lb > lc + ld) continue;
          bool ssss = false;
          if (la + lb + lc + ld == 0) ssss = true;
#if !SUPPORT_T1G12
          if (ssss) continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la, lb), max(lc, ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold =
              need_to_optimize ? cparams->unroll_threshold() : 0;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
          dg_xxxx->registry()->condense_expr(condense_expr(
              cparams->unroll_threshold(), cparams->max_vector_length() > 1));
          // Need to accumulate integrals?
          dg_xxxx->registry()->accumulate_targets(
              cparams->accumulate_targets());

          typedef R12kG12_11_11_sq int_type;
          typedef R12kG12 oper_type;
          // k=0
          if (!ssss) {
            typedef oper_type::Descriptor oper_descr;
            oper_type oper(oper_descr(0));
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            std::shared_ptr<int_type> abcd = int_type::Instance(
                *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // k=-1
          if (!ssss) {
            typedef oper_type::Descriptor oper_descr;
            oper_type oper(oper_descr(-1));
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            std::shared_ptr<int_type> abcd = int_type::Instance(
                *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

#if LIBINT_SUPPORT_T1G12
          // [T_1,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            typedef oper_type::Descriptor oper_descr;
            oper_type oper(oper_descr(0));
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            std::shared_ptr<int_type> abcd = int_type::Instance(
                *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          // [T_2,G12]
          if (true) {
            typedef TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            typedef oper_type::Descriptor oper_descr;
            oper_type oper(oper_descr(1));
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            std::shared_ptr<int_type> abcd = int_type::Instance(
                *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
#endif

          // [G12,[T1,G12]]
          if (!ssss) {
            typedef G12TiG12_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            typedef oper_type::Descriptor oper_descr;
            oper_type oper(
                oper_descr(0));  // doesn't matter whether T1 or T2 here
#if LIBINT_CONTRACTED_INTS
            oper.descr().contract();
#endif
            std::shared_ptr<int_type> abcd = int_type::Instance(
                *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          std::string _label;
          {
            ostringstream os;
            os << "(" << shells[la]->label() << " " << shells[lb]->label()
               << " | r_{12}^K * exp(-g*r_{12}^2) | " << shells[lc]->label()
               << " " << shells[ld]->label() << ")";
            _label = os.str();
          }

          std::string prefix(cparams->source_directory());
          std::string label(_label);
          std::deque<std::string> decl_filenames;
          std::deque<std::string> def_filenames;

          // this will generate code for this targets, and potentially
          // generate code for its prerequisites
          GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                       decl_filenames, def_filenames, prefix, label, false);

          // update max stack size
          const std::shared_ptr<TaskParameters>& tparams =
              taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(5);

          ostringstream oss;
          oss << context->label_to_function_name("libint2_build_r12kg12") << "["
              << la << "][" << lb << "][" << lc << "][" << ld
              << "] = " << context->label_to_function_name(label)
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          // need to declare this function internally
          for (std::deque<std::string>::const_iterator i =
                   decl_filenames.begin();
               i != decl_filenames.end(); ++i) {
            oss.str("");
            oss << "#include <" << *i << ">" << endl;
            iface->to_int_iface(oss.str());
          }

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxxx->reset();
          memman->reset();
        }  // end of d loop
      }    // end of c loop
    }      // end of b loop
  }        // end of a loop
}

#endif  // LIBINT_INCLUDE_G12

#ifdef LIBINT_INCLUDE_G12
void build_R12kG12_2b_2k_separate(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    std::shared_ptr<Libint2Iface>& iface) {
  // do not support this if the commutator integrals are needed
#if LIBINT_SUPPORT_T1G12
  assert(false);
#endif

  // Note that because r12_-1_g12 integrals are evaluated using the same RR as
  // ERI, no need to generate their code at all
  const int ntasks = 2;
  const char* task_names[] = {"r12_0_g12", "g12_T1_g12"};
  const char* task_NAMES[] = {"R12_0_R12", "G12_T1_G12"};

  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am("r12kg12");
  for (unsigned int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  for (int task = 0; task < ntasks; ++task) {
    const std::string task_name(task_names[task]);

    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    taskmgr.current(task_name);
    iface->to_params(
        iface->macro_define(std::string("MAX_AM_") + task_NAMES[task], lmax));
    iface->to_params(iface->macro_define("SUPPORT_T1G12", 0));

    std::shared_ptr<DirectedGraph> dg_xxxx(new DirectedGraph);
    std::shared_ptr<Strategy> strat(new Strategy);
    std::shared_ptr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
    std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
    std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());

    for (unsigned int la = 0; la <= lmax; la++) {
      for (unsigned int lb = 0; lb <= lmax; lb++) {
        for (unsigned int lc = 0; lc <= lmax; lc++) {
          for (unsigned int ld = 0; ld <= lmax; ld++) {
            if (la + lb + lc + ld == 0) continue;

            if (!ShellQuartetSetPredicate<static_cast<ShellSetType>(
                    LIBINT_SHELL_SET)>::value(la, lb, lc, ld))
              continue;

            using std::max;
            const unsigned int max_am = max(max(la, lb), max(lc, ld));
            const bool need_to_optimize =
                (max_am <= cparams->max_am_opt("r12kg12"));
            const unsigned int unroll_threshold =
                need_to_optimize ? cparams->unroll_threshold() : 0;
            dg_xxxx->registry()->unroll_threshold(unroll_threshold);
            dg_xxxx->registry()->do_cse(need_to_optimize);
            dg_xxxx->registry()->condense_expr(condense_expr(
                cparams->unroll_threshold(), cparams->max_vector_length() > 1));
            // Need to accumulate integrals?
            dg_xxxx->registry()->accumulate_targets(
                cparams->accumulate_targets());

            std::string _label;
            // k=0
            if (task == 0) {
              typedef R12kG12_11_11_sq int_type;
              typedef R12kG12 oper_type;
              typedef oper_type::Descriptor oper_descr;
              oper_type oper(oper_descr(0));
#if LIBINT_CONTRACTED_INTS
              oper.descr().contract();
#endif
              std::shared_ptr<int_type> abcd = int_type::Instance(
                  *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
              os << "building " << abcd->description() << endl;
              std::shared_ptr<DGVertex> abcd_ptr =
                  std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
              dg_xxxx->append_target(abcd_ptr);
              _label = abcd_ptr->label();
            }

            // [G12,[T1,G12]]
            if (task == 1) {
              typedef G12TiG12_11_11_sq int_type;
              typedef int_type::OperType oper_type;
              typedef oper_type::Descriptor oper_descr;
              oper_type oper(
                  oper_descr(0));  // doesn't matter whether T1 or T2 here
#if LIBINT_CONTRACTED_INTS
              oper.descr().contract();
#endif
              std::shared_ptr<int_type> abcd = int_type::Instance(
                  *shells[la], *shells[lb], *shells[lc], *shells[ld], 0u, oper);
              os << "building " << abcd->description() << endl;
              std::shared_ptr<DGVertex> abcd_ptr =
                  std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
              dg_xxxx->append_target(abcd_ptr);
              _label = abcd_ptr->label();
            }

            std::string prefix(cparams->source_directory());
            std::string label(_label);
            std::deque<std::string> decl_filenames;
            std::deque<std::string> def_filenames;

            // this will generate code for this targets, and potentially
            // generate code for its prerequisites
            GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman,
                         decl_filenames, def_filenames, prefix, label, false);

            // update max stack size
            const std::shared_ptr<TaskParameters>& tparams =
                taskmgr.current().params();
            tparams->max_stack_size(max_am, memman->max_memory_used());
            tparams->max_ntarget(1);

            ostringstream oss;
            oss << context->label_to_function_name(
                       std::string("libint2_build_") + task_names[task])
                << "[" << la << "][" << lb << "][" << lc << "][" << ld
                << "] = " << context->label_to_function_name(label)
                << context->end_of_stat() << endl;
            iface->to_static_init(oss.str());

            // need to declare this function internally
            for (std::deque<std::string>::const_iterator i =
                     decl_filenames.begin();
                 i != decl_filenames.end(); ++i) {
              oss.str("");
              oss << "#include <" << *i << ">" << endl;
              iface->to_int_iface(oss.str());
            }

#if DEBUG
            os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
            dg_xxxx->reset();
            memman->reset();
          }  // end of d loop
        }    // end of c loop
      }      // end of b loop
    }        // end of a loop
  }          // end of task loop
}

#endif  // LIBINT_INCLUDE_G12

#ifdef LIBINT_INCLUDE_G12DKH
void build_G12DKH_2b_2k(std::ostream& os,
                        const std::shared_ptr<CompilationParameters>& cparams,
                        std::shared_ptr<Libint2Iface>& iface) {
  const std::string task("g12dkh");
  vector<CGShell*> shells;
  unsigned int lmax = cparams->max_am(task);
  for (int l = 0; l <= lmax; l++) {
    shells.push_back(new CGShell(l));
  }
  ImplicitDimensions::set_default_dims(cparams);

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current(task);
  iface->to_params(iface->macro_define("MAX_AM_G12DKH", lmax));

  //
  // Construct graphs for each desired target integral and
  // 1) generate source code for the found traversal path
  // 2) extract all remaining unresolved recurrence relations and
  //    append them to the stack. Such unresolved RRs are RRs applied
  //    to sets of integrals (rather than to individual integrals).
  // 3) at the end, for each unresolved recurrence relation generate
  //    explicit source code
  //
  std::shared_ptr<DirectedGraph> dg_xxxx(new DirectedGraph);
  std::shared_ptr<Strategy> strat(new Strategy);
  std::shared_ptr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  for (int la = 0; la <= lmax; la++) {
    for (int lb = 0; lb <= lmax; lb++) {
      for (int lc = 0; lc <= lmax; lc++) {
        for (int ld = 0; ld <= lmax; ld++) {
          if (la < lb || lc < ld || la + lb > lc + ld) continue;
          bool ssss = false;
          if (la + lb + lc + ld == 0) ssss = true;

#if STUDY_MEMORY_USAGE
          const int lim = 5;
          if (!(la == lim && lb == lim && lc == lim && ld == lim)) continue;
#endif

          // unroll only if max_am <= cparams->max_am_opt(task)
          using std::max;
          const unsigned int max_am = max(max(la, lb), max(lc, ld));
          const bool need_to_optimize = (max_am <= cparams->max_am_opt(task));
          const unsigned int unroll_threshold =
              need_to_optimize ? cparams->unroll_threshold() : 0;
          dg_xxxx->registry()->unroll_threshold(unroll_threshold);
          dg_xxxx->registry()->do_cse(need_to_optimize);
          dg_xxxx->registry()->condense_expr(condense_expr(
              cparams->unroll_threshold(), cparams->max_vector_length() > 1));
          // Need to accumulate integrals?
          dg_xxxx->registry()->accumulate_targets(
              cparams->accumulate_targets());

          // k=0
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            typedef oper_type::Descriptor oper_descr;
            std::shared_ptr<int_type> abcd =
                int_type::Instance(*shells[la], *shells[lb], *shells[lc],
                                   *shells[ld], 0u, oper_type(oper_descr(0)));
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // k=2
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            typedef oper_type::Descriptor oper_descr;
            std::shared_ptr<int_type> abcd =
                int_type::Instance(*shells[la], *shells[lb], *shells[lc],
                                   *shells[ld], 0u, oper_type(oper_descr(2)));
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // k=4
          if (!ssss) {
            typedef R12kG12_11_11_sq int_type;
            typedef R12kG12 oper_type;
            typedef oper_type::Descriptor oper_descr;
            std::shared_ptr<int_type> abcd =
                int_type::Instance(*shells[la], *shells[lb], *shells[lc],
                                   *shells[ld], 0u, oper_type(oper_descr(4)));
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // (G12prime.Div1)^2
          if (true) {
            typedef DivG12prime_xTx_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            typedef oper_type::Descriptor oper_descr;
            std::shared_ptr<int_type> abcd =
                int_type::Instance(*shells[la], *shells[lb], *shells[lc],
                                   *shells[ld], 0u, oper_type(oper_descr(0)));
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }
          // (G12prime.Div2)^2
          if (true) {
            typedef DivG12prime_xTx_11_11_sq int_type;
            typedef int_type::OperType oper_type;
            typedef oper_type::Descriptor oper_descr;
            std::shared_ptr<int_type> abcd =
                int_type::Instance(*shells[la], *shells[lb], *shells[lc],
                                   *shells[ld], 0u, oper_type(oper_descr(1)));
            os << "building " << abcd->description() << endl;
            std::shared_ptr<DGVertex> abcd_ptr =
                std::dynamic_pointer_cast<DGVertex, int_type>(abcd);
            dg_xxxx->append_target(abcd_ptr);
          }

          std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
          std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());
          dg_xxxx->apply(strat, tactic);
          dg_xxxx->optimize_rr_out(context);
          dg_xxxx->traverse();
#if DEBUG
          os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
#endif

          std::string label;
          {
            ostringstream os;
            os << "(" << shells[la]->label() << " " << shells[lb]->label()
               << " | r_{12}^(0,2,4) * exp(-g*r_{12}^2) | "
               << shells[lc]->label() << " " << shells[ld]->label() << ")";
            label = os.str();
          }
          std::string prefix(cparams->source_directory());
          std::string decl_filename(prefix + context->label_to_name(label));
          decl_filename += ".h";
          std::string src_filename(prefix + context->label_to_name(label));
          src_filename += ".cc";
          std::basic_ofstream<char> declfile(decl_filename.c_str());
          std::basic_ofstream<char> srcfile(src_filename.c_str());
          dg_xxxx->generate_code(context, memman,
                                 ImplicitDimensions::default_dims(),
                                 std::shared_ptr<CodeSymbols>(new CodeSymbols),
                                 label, declfile, srcfile);

          // update max stack size
          const std::shared_ptr<TaskParameters>& tparams =
              taskmgr.current().params();
          tparams->max_stack_size(max_am, memman->max_memory_used());
          tparams->max_ntarget(3);

          ostringstream oss;
          oss << context->label_to_function_name("libint2_build_g12dkh") << "["
              << la << "][" << lb << "][" << lc << "][" << ld
              << "] = " << context->label_to_function_name(label)
              << context->end_of_stat() << endl;
          iface->to_static_init(oss.str());

          oss.str("");
          oss << "#include <" << decl_filename << ">" << endl;
          iface->to_int_iface(oss.str());

          // For the most expensive (i.e. presumably complete) graph extract
          // all precomputed quantities -- these will be members of the
          // evaluator structure also extract all RRs -- need to keep track of
          // these to figure out which external symbols appearing in RR code
          // belong to this task also
          if (la == lmax && lb == lmax && lc == lmax && ld == lmax)
            extract_symbols(dg_xxxx);

#if DEBUG
          os << "Max memory used = " << memman->max_memory_used() << endl;
#endif
          dg_xxxx->reset();
          declfile.close();
          srcfile.close();
        }  // end of d loop
      }    // end of c loop
    }      // end of b loop
  }        // end of a loop
}

#endif  // LIBINT_INCLUDE_G12DKH

void config_to_api(const std::shared_ptr<CompilationParameters>& cparams,
                   std::shared_ptr<Libint2Iface>& iface) {
  int max_deriv_order = 0;
#ifdef LIBINT_INCLUDE_ONEBODY
  iface->to_params(iface->macro_define("SUPPORT_ONEBODY", 1));
  iface->to_params(
      iface->macro_define("DERIV_ONEBODY_ORDER", LIBINT_INCLUDE_ONEBODY));
#ifdef LIBINT_DISABLE_ONEBODY_PROPERTY_DERIVS
  iface->to_params(iface->macro_define("DERIV_ONEBODY_PROPERTY_ORDER", 0));
#else
  iface->to_params(iface->macro_define("DERIV_ONEBODY_PROPERTY_ORDER",
                                       LIBINT_INCLUDE_ONEBODY));
#endif
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_ONEBODY);
#endif
#ifdef LIBINT_INCLUDE_ERI
  iface->to_params(iface->macro_define("SUPPORT_ERI", 1));
  iface->to_params(iface->macro_define("DERIV_ERI_ORDER", LIBINT_INCLUDE_ERI));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_ERI);
#endif
#ifdef LIBINT_INCLUDE_RKB_ERI
  iface->to_params(iface->macro_define("SUPPORT_RKB_ERI", 1));
  iface->to_params(
      iface->macro_define("DERIV_RKB_ERI_ORDER", LIBINT_INCLUDE_RKB_ERI));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_RKB_ERI);
#endif
#ifdef LIBINT_INCLUDE_ERI3
  iface->to_params(iface->macro_define("SUPPORT_ERI3", 1));
  iface->to_params(
      iface->macro_define("DERIV_ERI3_ORDER", LIBINT_INCLUDE_ERI3));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_ERI3);
#endif
#ifdef LIBINT_INCLUDE_ERI2
  iface->to_params(iface->macro_define("SUPPORT_ERI2", 1));
  iface->to_params(
      iface->macro_define("DERIV_ERI2_ORDER", LIBINT_INCLUDE_ERI2));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_ERI2);
#endif
#ifdef LIBINT_INCLUDE_G12
  iface->to_params(iface->macro_define("SUPPORT_G12", 1));
  iface->to_params(iface->macro_define("DERIV_G12_ORDER", LIBINT_INCLUDE_G12));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_G12);
#endif
#ifdef LIBINT_INCLUDE_G12DKH
  iface->to_params(iface->macro_define("SUPPORT_G12DKH", 1));
  iface->to_params(
      iface->macro_define("DERIV_G12DKH_ORDER", LIBINT_INCLUDE_G12DKH));
  max_deriv_order = std::max(max_deriv_order, LIBINT_INCLUDE_G12DKH);
#endif
  iface->to_params(iface->macro_define("MAX_DERIV_ORDER", max_deriv_order));

  // this is only needed for preprocessor-based generic processing of all
  // generated tasks declare all tasks in a range of valid tasks as defined or
  // not
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  // the range is defined by max # of centers, max deriv order, and operator
  // set
  const size_t max_ncenter = 4;
  for (unsigned int ncenter = 0; ncenter <= max_ncenter; ++ncenter) {
    std::stringstream oss;
    oss << ncenter;

    for (unsigned int d = 0; d <= max_deriv_order; ++d) {
      std::string abbrv_label, full_label;

      {  // 1-body ints
        std::string ncenter_str = oss.str();
        std::string ncenter_str_abbrv =
            ncenter == 2 ? std::string("") : oss.str();
#define BOOST_PP_MCR1(r, data, elem)                                         \
  abbrv_label = task_label(ncenter_str_abbrv + BOOST_PP_STRINGIZE(elem), d); \
  full_label = task_label(ncenter_str + BOOST_PP_STRINGIZE(elem), d);        \
  iface->to_params(                                                          \
      iface->macro_define(std::string("TASK_EXISTS_") + full_label,          \
                          taskmgr.exists(abbrv_label) ? 1 : 0));

        BOOST_PP_LIST_FOR_EACH(BOOST_PP_MCR1, _, BOOST_PP_ONEBODY_TASK_LIST)
#undef BOOST_PP_MCR1
      }

      {  // 2-body ints

#define BOOST_PP_TWOBODY_TASKOPER_TUPLE                                    \
  ("eri", "coulomb_opop", "opop_coulomb_opop", "op_coulomb_op", "r12kg12", \
   "r12_0_g12", "r12_2_g12", "g12_T1_g12", "g12dkh")
#define BOOST_PP_TWOBODY_TASKOPER_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_TWOBODY_TASKOPER_TUPLE)

        std::string ncenter_str = oss.str();
        std::string ncenter_str_abbrv =
            ncenter == 4 ? std::string("") : oss.str();
#define BOOST_PP_MCR1(r, data, elem)                                \
  abbrv_label = task_label(ncenter_str_abbrv + elem, d);            \
  full_label = task_label(ncenter_str + elem, d);                   \
  iface->to_params(                                                 \
      iface->macro_define(std::string("TASK_EXISTS_") + full_label, \
                          taskmgr.exists(abbrv_label) ? 1 : 0));

        BOOST_PP_LIST_FOR_EACH(BOOST_PP_MCR1, _, BOOST_PP_TWOBODY_TASKOPER_LIST)
#undef BOOST_PP_MCR1
      }
    }
  }
}
