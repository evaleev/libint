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

#include <algebra.h>
#include <code.h>
#include <context.h>
#include <dg.h>
#include <extract.h>
#include <graph_registry.h>
#include <integral.h>
#include <prefactors.h>
#include <rr.h>
#include <singl_stack.h>
#include <strategy.h>
#include <task.h>

#include <fstream>
#include <limits>

using namespace std;
using namespace libint2;
using namespace libint2::prefactor;

RecurrenceRelation::RecurrenceRelation() : nflops_(0), expr_() {}

RecurrenceRelation::~RecurrenceRelation() {}

//
// If there is no generic equivalent, generate explicit code for this recurrence
// relation: 1) append target and children to a DirectedGraph dg 2) set their
// code symbols 3) apply IntSet_to_Ints 4) Apply RRs such that no additional
// vertices appear 5) call dg->generate_code()
//
void RecurrenceRelation::generate_code(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims,
    const std::string& funcname, std::ostream& decl, std::ostream& def) {
  //
  // Check if there is a generic equivalent that can be used
  //
  if (this->has_generic(context->cparams())) {
    generate_generic_code(context, dims, funcname, decl, def);
    return;
  }

  const std::shared_ptr<DGVertex> target_vptr = rr_target();
#if DEBUG
  std::cout << "RecurrenceRelation::generate_code: target = "
            << target_vptr->label() << std::endl;
#endif

  const std::shared_ptr<CompilationParameters>& cparams = context->cparams();
  std::shared_ptr<DirectedGraph> dg(new DirectedGraph);
  dg->set_label(context->label_to_name(label_to_funcname(funcname)));
  generate_graph_(dg);

  // Intermediates in RR code are either are automatic variables or have to go
  // on vstack
  dg->registry()->stack_name("inteval->vstack");
  // No need to return the targets via inteval's targets
  dg->registry()->return_targets(false);

  // check if CSE to be performed
  typedef IntegralSet<IncableBFSet> ISet;
  std::shared_ptr<ISet> target =
      std::dynamic_pointer_cast<ISet, DGVertex>(target_vptr);
  if (target) {
    //
    // do CSE only if max_am <= cparams->max_am_opt()
    //
    const unsigned int np = target->num_part();
    unsigned int max_am = 0;
    // bra
    for (unsigned int p = 0; p < np; p++) {
      const unsigned int nf = target->num_func_bra(p);
      for (unsigned int f = 0; f < nf; f++) {
        // Assuming shells here
        const unsigned int am = target->bra(p, f).norm();
        using std::max;
        max_am = max(max_am, am);
      }
    }
    // ket
    for (unsigned int p = 0; p < np; p++) {
      const unsigned int nf = target->num_func_ket(p);
      for (unsigned int f = 0; f < nf; f++) {
        // Assuming shells here
        const unsigned int am = target->ket(p, f).norm();
        using std::max;
        max_am = max(max_am, am);
      }
    }
    const bool need_to_optimize = (max_am <= cparams->max_am_opt());
    dg->registry()->do_cse(need_to_optimize);
  }
  dg->registry()->condense_expr(
      condense_expr(std::numeric_limits<unsigned int>::max(),
                    cparams->max_vector_length() > 1));
  dg->registry()->ignore_missing_prereqs(
      true);  // assume all prerequisites are available -- if some are not,
              // something is VERY broken

#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".strat.dot");
    dg->print_to_dot(false, dotfile);
  }
#endif

  // Assign symbols for the target and source integral sets
  std::shared_ptr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);
  // Traverse the graph
  dg->optimize_rr_out(context);
  dg->traverse();
#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".expr.dot");
    dg->print_to_dot(false, dotfile);
  }
#endif
  // Generate code
  std::shared_ptr<MemoryManager> memman(new WorstFitMemoryManager());
  std::shared_ptr<ImplicitDimensions> localdims = adapt_dims_(dims);
  dg->generate_code(context, memman, localdims, symbols, funcname, decl, def);

  // extract all external symbols -- these will be members of the evaluator
  // structure
  std::shared_ptr<ExtractExternSymbols> extractor(new ExtractExternSymbols);
  dg->foreach (*extractor);
  const ExtractExternSymbols::Symbols& externsymbols = extractor->symbols();

#if 0
  // print out the symbols
  std::cout << "Recovered symbols from DirectedGraph for " << label() << std::endl;
  typedef ExtractExternSymbols::Symbols::const_iterator citer;
  citer end = externsymbols.end();
  for(citer t=externsymbols.begin(); t!=end; ++t)
    std::cout << *t << std::endl;
#endif

#if PRINT_DAG_GRAPHVIZ
  {
    std::basic_ofstream<char> dotfile(dg->label() + ".symb.dot");
    dg->print_to_dot(false, dotfile);
  }
#endif

  // get this RR InstanceID
  RRStack::InstanceID myid =
      RRStack::Instance()
          ->find(std::enable_shared_from_this<this_type>::shared_from_this())
          .first;

  // For each task which requires this RR:
  // 1) update max stack size
  // 2) append external symbols from this RR to its list
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  for (tciter t = taskmgr.first(); t != tend; ++t) {
    const std::shared_ptr<TaskExternSymbols> tsymbols = t->symbols();
    if (tsymbols->find(myid)) {
      // update max stack size
      t->params()->max_vector_stack_size(memman->max_memory_used());
      // add external symbols
      tsymbols->add(externsymbols);
    }
  }

  dg->reset();
}

namespace libint2 {
// generate_generic_code reuses this function from dg.cc:
extern std::string declare_function(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims,
    const std::shared_ptr<CodeSymbols>& args, const std::string& tlabel,
    const std::string& function_descr, std::ostream& decl);
}  // namespace libint2

void RecurrenceRelation::generate_generic_code(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims,
    const std::string& funcname, std::ostream& decl, std::ostream& def) {
  const std::shared_ptr<DGVertex> target_vptr = rr_target();
  std::cout << "RecurrenceRelation::generate_generic_code: target = "
            << target_vptr->label() << std::endl;

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  const std::string tlabel = taskmgr.current().label();
  std::shared_ptr<ImplicitDimensions> localdims = adapt_dims_(dims);
  // Assign symbols for the target and source integral sets
  std::shared_ptr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);

  // declare function
  const std::string func_decl =
      declare_function(context, localdims, symbols, tlabel, funcname, decl);

  //
  // Generate function's definition
  //

  // include standard headers
  def << context->std_header();
  //         + generic code declaration
  def << "#include <" << this->generic_header() << ">" << endl;
  def << endl;

  // start the body ...
  def << context->code_prefix();
  def << func_decl << context->open_block() << endl;
  def << context->std_function_header();

  // ... fill the body
  def << this->generic_instance(context, symbols) << endl;

  // ... end the body
  def << context->close_block() << endl;
  def << context->code_postfix();
}

std::shared_ptr<DirectedGraph> RecurrenceRelation::generate_graph_(
    const std::shared_ptr<DirectedGraph>& dg) {
  dg->append_target(rr_target());
  for (unsigned int c = 0; c < num_children(); c++)
    dg->append_vertex(rr_child(c));
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets = "
       << dg->num_vertices() << endl;
#endif
  std::shared_ptr<Strategy> strat(new Strategy);
  std::shared_ptr<Tactic> ntactic(new NullTactic);
  // Always need to unroll integral sets first
  dg->registry()->unroll_threshold(std::numeric_limits<unsigned int>::max());
  dg->apply(strat, ntactic);
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets + "
          "integrals = "
       << dg->num_vertices() << endl;
#endif
  // Mark children sets and their descendants to not compute
  for (unsigned int c = 0; c < num_children(); c++)
    dg->apply_at<&DGVertex::not_need_to_compute>(rr_child(c));
  // Apply recurrence relations using existing vertices on the graph (i.e.
  // such that no new vertices appear)
  std::shared_ptr<Tactic> ztactic(new FewestNewVerticesTactic(dg));
  dg->apply(strat, ztactic);
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- should be same as previous = "
       << dg->num_vertices() << endl;
#endif

  return dg;
}

void RecurrenceRelation::assign_symbols_(
    std::shared_ptr<CodeSymbols>& symbols) {
  // Set symbols on the target and children sets
  rr_target()->set_symbol("target");
  symbols->append_symbol("target");
  for (unsigned int c = 0; c < num_children(); c++) {
    ostringstream oss;
    oss << "src" << c;
    string symb = oss.str();
    rr_child(c)->set_symbol(symb);
    symbols->append_symbol(symb);
  }
}

std::shared_ptr<ImplicitDimensions> RecurrenceRelation::adapt_dims_(
    const std::shared_ptr<ImplicitDimensions>& dims) const {
  return dims;
}

const std::string& RecurrenceRelation::label() const {
  if (label_.empty()) label_ = generate_label();
  return label_;
}

std::string RecurrenceRelation::description() const {
  const std::string descr = label();
  return descr;
}

void RecurrenceRelation::add_expr(const std::shared_ptr<ExprType>& expr,
                                  int minus) {
  if (expr_ == 0) {
    if (minus != -1) {
      expr_ = expr;
    } else {
      using libint2::prefactor::Scalar;
      std::shared_ptr<ExprType> negative(
          new ExprType(ExprType::OperatorTypes::Times, expr, Scalar(-1.0)));
      expr_ = negative;
      ++nflops_;
    }
  } else {
    if (minus != -1) {
      std::shared_ptr<ExprType> sum(
          new ExprType(ExprType::OperatorTypes::Plus, expr_, expr));
      expr_ = sum;
      ++nflops_;
    } else {
      std::shared_ptr<ExprType> sum(
          new ExprType(ExprType::OperatorTypes::Minus, expr_, expr));
      expr_ = sum;
      ++nflops_;
    }
  }
}

bool RecurrenceRelation::invariant_type() const {
  // By default, recurrence relations do not change the type of the functions,
  // i.e. VRR applied to an integral over shells will produce integrals over
  // shells
  return true;
}

std::string RecurrenceRelation::spfunction_call(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims) const {
  ostringstream os;
  os << context->label_to_name(
            label_to_funcname(context->cparams()->api_prefix() + label()))
     // First argument is the library object
     << "(inteval, "
     // Second is the target
     << context->value_to_pointer(rr_target()->symbol());
  // then come children
  const unsigned int nchildren = num_children();
  for (unsigned int c = 0; c < nchildren; c++) {
    os << ", " << context->value_to_pointer(rr_child(c)->symbol());
  }
  os << ")" << context->end_of_stat() << endl;
  return os.str();
}

bool RecurrenceRelation::has_generic(
    const std::shared_ptr<CompilationParameters>& cparams) const {
  return false;
}

std::string RecurrenceRelation::generic_header() const {
  throw std::logic_error(
      "RecurrenceRelation::generic_header() -- should not be called! Check if "
      "DerivedRecurrenceRelation::generic_header() is implemented");
}

std::string RecurrenceRelation::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  throw std::logic_error(
      "RecurrenceRelation::generic_instance() -- should not be called! Check "
      "if DerivedRecurrenceRelation::generic_instance() is implemented");
}

size_t RecurrenceRelation::size_of_children() const {
  const auto nchildren = this->num_children();
  size_t result = 0;
  for (auto c = 0; c != nchildren; ++c) {
    result += this->rr_child(c)->size();
  }
  return result;
}

namespace libint2 {
namespace algebra {
/// these operators are extremely useful to write compact expressions
std::shared_ptr<RecurrenceRelation::ExprType> operator+(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  return std::shared_ptr<Oper>(new Oper(Oper::OperatorTypes::Plus, A, B));
}
std::shared_ptr<RecurrenceRelation::ExprType> operator-(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  return std::shared_ptr<Oper>(new Oper(Oper::OperatorTypes::Minus, A, B));
}
std::shared_ptr<RecurrenceRelation::ExprType> operator*(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  return std::shared_ptr<Oper>(new Oper(Oper::OperatorTypes::Times, A, B));
}
std::shared_ptr<RecurrenceRelation::ExprType> operator/(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  return std::shared_ptr<Oper>(new Oper(Oper::OperatorTypes::Divide, A, B));
}
const std::shared_ptr<RecurrenceRelation::ExprType>& operator+=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  if (A) {
    const std::shared_ptr<Oper>& Sum = A + B;
    A = Sum;
  } else
    A = Scalar(0) + B;
  return A;
}
const std::shared_ptr<RecurrenceRelation::ExprType>& operator-=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  if (A) {
    const std::shared_ptr<Oper>& Diff = A - B;
    A = Diff;
  } else
    A = Scalar(0) - B;
  return A;
}
const std::shared_ptr<RecurrenceRelation::ExprType>& operator*=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  const std::shared_ptr<Oper>& Product = A * B;
  A = Product;
  return A;
}
const std::shared_ptr<RecurrenceRelation::ExprType>& operator/=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B) {
  typedef RecurrenceRelation::ExprType Oper;
  const std::shared_ptr<Oper>& Quotient = A / B;
  A = Quotient;
  return A;
}
}  // namespace algebra
}  // namespace libint2

///////////////
