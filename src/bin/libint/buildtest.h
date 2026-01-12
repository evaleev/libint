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

#ifndef _libint2_src_bin_libint_buildtest_h_
#define _libint2_src_bin_libint_buildtest_h_

#include <dg.h>
#include <dims.h>
#include <graph_registry.h>
#include <iface.h>
#include <integral_11_11.h>
#include <strategy.h>

#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

namespace libint2 {

// defined in buildtest.cc
void generate_rr_code(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams,
                      std::deque<std::string>& decl_filenames,
                      std::deque<std::string>& def_filenames);

/// defined below generates code for dg; dg and memman are reset at the end
void GenerateCode(const std::shared_ptr<DirectedGraph>& dg,
                  const std::shared_ptr<CodeContext>& context,
                  const std::shared_ptr<CompilationParameters>& cparams,
                  const std::shared_ptr<Strategy>& strat,
                  const std::shared_ptr<Tactic>& tactic,
                  const std::shared_ptr<MemoryManager>& memman,
                  std::deque<std::string>& decl_filenames,
                  std::deque<std::string>& def_filenames,
                  const std::string& prefix, const std::string& label,
                  bool have_parent);

/// Command-line parser for the standard build tester -- N is the number of
/// centers, i.e. 4 for 4-center ERI
template <unsigned int N>
class TesterCmdLine {
 public:
  TesterCmdLine(int argc, char* argv[]);
  ~TesterCmdLine() {}

  const std::vector<unsigned int>& am() const { return am_; }
  unsigned int size_to_unroll() const { return size_to_unroll_; }
  unsigned int veclen() const { return veclen_; }
  bool vectorize_by_line() const { return vectorize_by_line_; }
  bool do_cse() const { return do_cse_; }

 private:
  static const unsigned int max_am = 10;
  std::vector<unsigned int> am_;
  unsigned int size_to_unroll_;
  unsigned int veclen_;
  bool vectorize_by_line_;
  bool do_cse_;
};

/** This is a user-friendly generic test of building an Integral using specified
   size_to_unroll, veclen, vec_by_line, and do_cse. GenAllCode should be set to
   true if compilable code to be produced (i.e. include header files + set-level
   recurrence relations code)
 */
template <class Integral, bool GenAllCode>
void BuildTest(const std::vector<std::shared_ptr<Integral> >& targets,
               unsigned int size_to_unroll, unsigned int veclen,
               bool vec_by_line, bool do_cse,
               const std::string& complabel = "buildtest",
               std::ostream& os = std::cout);

/** This is a generic test of building an Integral using specified cparams,
   memman, size_to_unroll, default strategy and specified tactic. GenAllCode
   should be set to true if compilable code to be produced (i.e. include header
   files + set-level recurrence relations code)
 */
template <class Integral, bool GenAllCode>
void __BuildTest(
    const std::vector<std::shared_ptr<Integral> >& targets,
    const std::shared_ptr<CompilationParameters>& cparams,
    unsigned int size_to_unroll, std::ostream& os = std::cout,
    const std::shared_ptr<Tactic>& tactic =
        std::shared_ptr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
    const std::shared_ptr<MemoryManager>& memman =
        std::shared_ptr<MemoryManager>(new WorstFitMemoryManager),
    const std::string& complabel = "general_integral");

template <class Integral, bool GenAllCode>
void __BuildTest(const std::vector<std::shared_ptr<Integral> >& targets,
                 const std::shared_ptr<CompilationParameters>& cparams,
                 unsigned int size_to_unroll, std::ostream& os,
                 const std::shared_ptr<Tactic>& tactic,
                 const std::shared_ptr<MemoryManager>& memman,
                 const std::string& complabel) {
  const std::string prefix("");
  const std::string label = cparams->api_prefix() + complabel;
  std::shared_ptr<Strategy> strat(new Strategy);
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.add(complabel);
  taskmgr.current(complabel);

  //
  // do CSE only if max_am <= cparams->max_am_opt()
  //
  unsigned int max_am = 0;
  for (unsigned int t = 0; t < targets.size(); ++t) {
    const std::shared_ptr<Integral>& target = targets[t];
    const unsigned int np = target->bra().num_part();
    // bra
    for (unsigned int p = 0; p < np; p++) {
      const unsigned int nf = target->bra().num_members(p);
      for (unsigned int f = 0; f < nf; f++) {
        // Assuming shells here
        const unsigned int am = target->bra(p, f).qn();
        using std::max;
        max_am = max(max_am, am);
      }
    }
    // ket
    for (unsigned int p = 0; p < np; p++) {
      const unsigned int nf = target->ket().num_members(p);
      for (unsigned int f = 0; f < nf; f++) {
        // Assuming shells here
        const unsigned int am = target->ket(p, f).qn();
        using std::max;
        max_am = max(max_am, am);
      }
    }
  }
  const bool need_to_optimize = (max_am <= cparams->max_am_opt(complabel));

  std::deque<std::string> decl_filenames;
  std::deque<std::string> def_filenames;

  os << "Building " << complabel << std::endl;

  std::shared_ptr<DirectedGraph> dg_xxxx(new DirectedGraph);
  dg_xxxx->set_label(complabel);

  // configure the graph
  dg_xxxx->registry()->do_cse(need_to_optimize);
  dg_xxxx->registry()->condense_expr(
      condense_expr(size_to_unroll, cparams->max_vector_length() > 1));
  // Need to accumulate integrals?
  dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());
  dg_xxxx->registry()->unroll_threshold(size_to_unroll);

  for (unsigned int t = 0; t < targets.size(); ++t) {
    const std::shared_ptr<Integral>& target = targets[t];
    std::shared_ptr<DGVertex> target_ptr =
        std::dynamic_pointer_cast<DGVertex, Integral>(target);
    assert(target_ptr != 0);
    dg_xxxx->append_target(target_ptr);
  }

  // this will generate code for this targets, and potentially generate code for
  // its prerequisites
  GenerateCode(dg_xxxx, context, cparams, strat, tactic, memman, decl_filenames,
               def_filenames, prefix, label, false);

  // update max stack size
  taskmgr.current().params()->max_stack_size(max_am, memman->max_memory_used());
  taskmgr.current().params()->max_ntarget(targets.size());
  os << "Max memory used = " << memman->max_memory_used() << std::endl;

  if (GenAllCode) {
    // initialize code context to produce library API
    std::shared_ptr<CodeContext> icontext(new CppCodeContext(cparams));
    // initialize object to generate interface
    std::shared_ptr<Libint2Iface> iface(new Libint2Iface(cparams, icontext));

    // generate interface
    std::ostringstream oss;
    for (std::deque<std::string>::const_iterator i = decl_filenames.begin();
         i != decl_filenames.end(); ++i) {
      oss << "#include <" << *i << ">" << std::endl;
    }
    iface->to_int_iface(oss.str());

    // transfer some configuration parameters to the generated library API
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

    // Generate set-level RR code
    generate_rr_code(os, cparams, decl_filenames, def_filenames);

    // Print log
    std::cout << "Generated headers: ";
    std::copy(decl_filenames.begin(), decl_filenames.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl << "Generated sources: ";
    std::copy(def_filenames.begin(), def_filenames.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl
              << "Top compute function: "
              << context->label_to_name(label_to_funcname(label)) << std::endl;
  }
}

void GenerateCode(const std::shared_ptr<DirectedGraph>& dg,
                  const std::shared_ptr<CodeContext>& context,
                  const std::shared_ptr<CompilationParameters>& cparams,
                  const std::shared_ptr<Strategy>& strat,
                  const std::shared_ptr<Tactic>& tactic,
                  const std::shared_ptr<MemoryManager>& memman,
                  std::deque<std::string>& decl_filenames,
                  std::deque<std::string>& def_filenames,
                  const std::string& prefix, const std::string& label,
                  bool have_parent) {
  dg->apply(strat, tactic);
#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".strat.dot");
    dg->print_to_dot(false, dotfile);
  }
#endif
  dg->optimize_rr_out(context);
#if DEBUG
  std::cout << "The number of vertices = " << dg->num_vertices() << std::endl;
#endif

  // if there are missing prerequisites -- make a list of them
  PrerequisitesExtractor pe;
  if (dg->missing_prerequisites()) {
    // std::cout << "missing some prerequisites!" << std::endl;
    dg->foreach (pe);
  }
  std::deque<std::shared_ptr<DGVertex> > prereq_list = pe.vertices;

  dg->traverse();
  // dg->debug_print_traversal(cout);

#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".expr.dot");
    dg->print_to_dot(false, dotfile);
  }
#endif

  std::string decl_filename(prefix + context->label_to_name(label));
  decl_filename += ".h";
  std::string def_filename(prefix + context->label_to_name(label));
  def_filename += ".cc";
  std::basic_ofstream<char> declfile(decl_filename.c_str());
  std::basic_ofstream<char> deffile(def_filename.c_str());
  // if have parent graph, it will pass its stack where this graph will put its
  // results
  std::shared_ptr<CodeSymbols> args(new CodeSymbols);
  if (have_parent) args->append_symbol("parent_stack");
  dg->generate_code(context, memman, ImplicitDimensions::default_dims(), args,
                    label, declfile, deffile);
  declfile.close();
  deffile.close();

  // extract all external symbols
  extract_symbols(dg);

#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".symb.dot");
    dg->print_to_dot(true, dotfile);
  }
#endif

  decl_filenames.push_back(decl_filename);
  def_filenames.push_back(def_filename);

  // last: missing prerequisites? create new graph computing prereqs and move
  // them onto it
  if (dg->missing_prerequisites()) {
    std::shared_ptr<DirectedGraph> dg_prereq(new DirectedGraph);
    // configure identically
    dg_prereq->registry() =
        std::shared_ptr<GraphRegistry>(dg->registry()->clone());
    // except:
    // - allow uncontraction
    // - no need to return targets via inteval->targets_ -- their locations are
    // known by the parent graph (see allocate_mem)
    dg_prereq->registry()->uncontract(true);
    assert(cparams->contracted_targets());
    dg_prereq->registry()->return_targets(false);
    dg_prereq->registry()->accumulate_targets(true);
    dg_prereq->registry()->stack_name("stack");
    if (dg->registry()->current_timer() >= 0) {
      dg_prereq->registry()->current_timer(dg->registry()->current_timer() + 1);
    }

    // now is the "right" time to reset dg
    // reset graph of the previous computation so that the vertices that will be
    // targets on the new graph are not attached still to the vertices from the
    // old graph
    dg->reset();
    memman->reset();

    while (!prereq_list.empty()) {
      dg_prereq->append_target(prereq_list.front());
      prereq_list.pop_front();
    }

    const std::string label_prereq = label + "_prereq";
    GenerateCode(dg_prereq, context, cparams, strat, tactic, memman,
                 decl_filenames, def_filenames, prefix, label_prereq, true);
  }
  dg->reset();
  memman->reset();
}

template <class Integral, bool GenAllCode>
void BuildTest(const std::vector<std::shared_ptr<Integral> >& targets,
               unsigned int size_to_unroll, unsigned int veclen,
               bool vec_by_line, bool do_cse, const std::string& complabel,
               std::ostream& os) {
  const unsigned int max_am = 10;
  os << "generating code to compute " << complabel << std::endl;

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.add(complabel);
  taskmgr.current(complabel);

  // initialize cparams
  std::shared_ptr<CompilationParameters> cparams(new CompilationParameters);
  cparams->max_am(complabel, max_am);
  cparams->num_bf(complabel, 4u);
  cparams->max_vector_length(veclen);
  cparams->vectorize_by_line(vec_by_line);
#if LIBINT_ALIGN_SIZE
  cparams->align_size(LIBINT_ALIGN_SIZE);
#endif
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
#if LIBINT_CONTRACTED_INTS
  cparams->contracted_targets(true);
#else
  cparams->contracted_targets(false);
#endif
#ifdef LIBINT_USER_DEFINED_REAL
  {
    const std::string realtype(LIBINT_USER_DEFINED_REAL);
    cparams->realtype(realtype);
  }
#endif

  if (do_cse) {
    cparams->max_am_opt(complabel, max_am);
  } else {
    cparams->max_am_opt(complabel, 0);
  }
  cparams->default_task_name(complabel);

  // set default dims
  ImplicitDimensions::set_default_dims(cparams);

  std::shared_ptr<StdRandomizePolicy> rpolicy(new StdRandomizePolicy(0.00));
  // use 4-center OS if the target is a 4-center integral
  std::shared_ptr<Tactic> tactic;
  {
    typedef GenIntegralSet_11_11<typename Integral::BasisFunctionType,
                                 typename Integral::OperatorType,
                                 typename Integral::AuxIndexType>
        genint_11_11_t;
    std::shared_ptr<genint_11_11_t> cast_ptr =
        std::dynamic_pointer_cast<genint_11_11_t>(targets.front());
    if (cast_ptr) {
      const unsigned int la = cast_ptr->bra(0, 0).norm();
      const unsigned int lb = cast_ptr->ket(0, 0).norm();
      const unsigned int lc = cast_ptr->bra(1, 0).norm();
      const unsigned int ld = cast_ptr->ket(1, 0).norm();
      tactic =
          std::shared_ptr<Tactic>(new FourCenter_OS_Tactic(la, lb, lc, ld));
    } else {
      tactic = std::shared_ptr<Tactic>(
          new FirstChoiceTactic<StdRandomizePolicy>(rpolicy));
    }
  }
  const std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager);
  __BuildTest<Integral, true>(targets, cparams, size_to_unroll, os, tactic,
                              memman, complabel);
}

template <unsigned int N>
TesterCmdLine<N>::TesterCmdLine(int argc, char* argv[]) {
  if (N == 0)
    throw ProgrammingError("TesterCmdLine<N>::TesterCmdLine but N is 0");
  const int argc_min = N + 2;
  const int argc_max = N + 5;
  if (argc < argc_min || argc > argc_max) {
    std::cerr
        << "Usage: " << argv[0]
        << " <am> size_to_unroll [vector_length] [vector_method] [do_cse]"
        << std::endl
        << "       <am> -- angular momenta on each center, e.g. 4 nonnegative "
           "integers for a 4-center ERI"
        << std::endl
        << "       size_to_unroll -- size of the largest integral set to be "
           "unrolled"
        << std::endl
        << "       vector_length  -- (optional) max vector length. Defaults to "
           "1."
        << std::endl
        << "       vector_method  -- (optional) vectorization method. Valid "
           "choices are 0 (by-block) and 1 (by-line). Defaults to 0."
        << std::endl
        << "       do_cse  -- (optional) do Common Subexpression Elimination? "
           "Valid choices are 0 (no) and 1 (yes). Defaults to 0."
        << std::endl
        << std::endl;
    throw InputError(
        "TesterCmdLine<N>::TesterCmdLine -- incorrect number of command-line "
        "arguments");
  }
  for (unsigned int i = 1; i < N + 1; ++i) {
    const unsigned int am = atoi(argv[i]);
    if (am > max_am)
      throw InputError(
          "TesterCmdLine<N>::TesterCmdLine -- angular momentum limit exceeded");
    am_.push_back(am);
  }
  size_to_unroll_ = atoi(argv[N + 1]);

  veclen_ = 1;
  if (argc >= N + 3) {
    veclen_ = atoi(argv[N + 2]);
  }
  vectorize_by_line_ = false;
  if (argc >= N + 4) {
    vectorize_by_line_ = (1 == atoi(argv[N + 3]));
  }
  do_cse_ = false;
  if (argc >= N + 5) {
    do_cse_ = (1 == atoi(argv[N + 4]));
  }
}

};  // namespace libint2

#endif
