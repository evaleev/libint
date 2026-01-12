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

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

#include <bfset.h>
#include <code.h>
#include <default_params.h>
#include <exception.h>
#include <global_macros.h>
#include <polyconstr.h>
#include <singl_stack.h>
#include <smart_ptr.h>
#include <util_types.h>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

class DGVertex;
class CodeContext;
class ImplicitDimensions;
class DirectedGraph;
template <typename V>
class AlgebraicOperator;

/** RRStack implements a stack of RecurrenceRelation's which can only hold
    one instance of a given RR. RecurrenceRelation::label() is used for hashing
  */
template <typename RR>
class RRStackBase : public SingletonStack<RR, std::string> {
 private:
  static std::shared_ptr<RRStackBase<RR> > rrstack_;

  // private constructor because it's a Singleton
  RRStackBase() : parent_type(&RR::label) {}

 public:
  typedef SingletonStack<RR, std::string> parent_type;
  typedef typename parent_type::data_type data_type;
  typedef typename parent_type::value_type value_type;
  typedef typename parent_type::iter_type iter_type;
  typedef typename parent_type::citer_type citer_type;
  typedef typename parent_type::InstanceID InstanceID;

  /// Obtain the unique Instance of RRStack
  static std::shared_ptr<RRStackBase<RR> >& Instance() {
    if (!rrstack_) {
      std::shared_ptr<RRStackBase<RR> > tmpstack(new RRStackBase<RR>);
      rrstack_ = tmpstack;
    }
    return rrstack_;
  }

  virtual ~RRStackBase() {}

  /// adds content of rrs to this stack
  void add(const std::shared_ptr<RRStackBase<RR> >& rrs) {
    for (citer_type it = rrs->begin(); it != rrs->end(); it++) {
      find((*it).second.second);
    }
  }

  /// removes rr from the stack
  void remove(const data_type& rr) { parent_type::remove(rr); }
};

template <typename RR>
std::shared_ptr<RRStackBase<RR> > RRStackBase<RR>::rrstack_;

/**
   RecurrenceRelation describes all recurrence relations
*/
class RecurrenceRelation
    : public std::enable_shared_from_this<RecurrenceRelation> {
 public:
  typedef RecurrenceRelation this_type;

  RecurrenceRelation();
  virtual ~RecurrenceRelation();

  /** num_children() returns the actual number of children.
      For example, VRR for ERIs has 5 children on the right-hand side,
      however, for some ERI classes (ss|ps) the actual number may be
      smaller.
  */
  virtual unsigned int num_children() const = 0;
  /// Returns i-th child
  virtual std::shared_ptr<DGVertex> rr_child(unsigned int i) const = 0;
  /// Returns the target
  virtual std::shared_ptr<DGVertex> rr_target() const = 0;

  /** Numerical expression of a recurrence relation is always expressed as
      an AlgebraicOperator<DGVertex> */
  typedef AlgebraicOperator<DGVertex> ExprType;
  /// Returns the expression
  const std::shared_ptr<ExprType>& rr_expr() const { return expr_; }

  /**
     Returns true is this recurrence relation is simple enough to optimize away.
     As a result of such optimization, standalone function will NOT be
     generated for this recurrence relation. Instead, it's source will be
     inlined and optimized.
  */
  virtual bool is_simple() const = 0;
  /**
     Returns true is the type of target and all children are exactly the same
  */
  virtual bool invariant_type() const;
  /**
   * @return 1 if recurrence relation transfers quanta from particle \c from to
   * particle \c to  where \c from < \c to , -1 if \c from > \c to , and 0 if
   * neither
   */
  virtual int partindex_direction() const { return 0; }
  /**
   * @return BraketDirection::BraToKet if recurrence relation transfers quanta
   * from function in bra to function in ket, BraketDirection::KetToBra if the
   * transfer is from ket to bra, and BraketDirection::None if neither.
   */
  virtual BraketDirection braket_direction() const {
    return BraketDirection::None;
  }
  /**
   * @return the total size of the children of this RR
   */
  size_t size_of_children() const;

  /**
    label() returns a unique, short, descriptive label of this RR
    (e.g. "VRR A (p s | 1/r_{12} | d s )" for Obara-Saika recurrence relation
    applied to center A to compute (ps|ds) ERI)
  */
  const std::string& label() const;

  /**
    description() returns a verbose description of this RR
  */
  virtual std::string description() const;

  /// Generate declaration and definition for the recurrence relation
  virtual void generate_code(const std::shared_ptr<CodeContext>& context,
                             const std::shared_ptr<ImplicitDimensions>& dims,
                             const std::string& funcname, std::ostream& decl,
                             std::ostream& def);

  /// Generate declaration and definition for the recurrence relation
  /// using generic code (typically, a manually written code)
  virtual void generate_generic_code(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims,
      const std::string& funcname, std::ostream& decl, std::ostream& def);

  /// Generate a callback for this recurrence relation
  virtual std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const;

  /// Return the number of FLOPs per this recurrence relation
  unsigned int nflops() const { return nflops_; }

  /// RecurrenceRelation is managed by SingletonStack but doesn't need to keep
  /// track of instance ID
  void inst_id(
      const SingletonStack<RecurrenceRelation, std::string>::InstanceID& i) {}

 protected:
  unsigned int nflops_;
  mutable std::string label_;
  std::shared_ptr<ExprType> expr_;
  /// Adds a (or -a, if minus = -1) to expr_.
  void add_expr(const std::shared_ptr<ExprType>& a, int minus = 1);
  /// Generates the label
  virtual std::string generate_label() const = 0;
  /// Registers with the stack
  template <class RR>
  bool register_with_rrstack() {
    // only register RRs with for shell sets
    if (TrivialBFSet<typename RR::BasisFunctionType>::result) return false;
    std::shared_ptr<RRStackBase<RecurrenceRelation> > rrstack =
        RRStackBase<RecurrenceRelation>::Instance();
    std::shared_ptr<RR> this_ptr = std::const_pointer_cast<RR, const RR>(
        std::static_pointer_cast<const RR, const RecurrenceRelation>(
            std::enable_shared_from_this<
                RecurrenceRelation>::shared_from_this()));
    rrstack->find(this_ptr);
#if DEBUG || DEBUG_CONSTRUCTION
    std::cout << "register_with_rrstack: registered " << this_ptr->label()
              << std::endl;
#endif
    return true;
  }

 private:
  /** used by generate_code to initialize the computation graph that computes
   * sets of integrals using the RR
   */
  std::shared_ptr<DirectedGraph> generate_graph_(
      const std::shared_ptr<DirectedGraph>& dg);
  /** assigns "target" symbol to the target vertex and "src<i>" to the i-th
     child vertex. Also appends these symbols to S. */
  void assign_symbols_(std::shared_ptr<CodeSymbols>& S);
  /** given an ImplicitDimension for the computation, adapt it for this
     recurrence relation. Default version does not do anything. */
  virtual std::shared_ptr<ImplicitDimensions> adapt_dims_(
      const std::shared_ptr<ImplicitDimensions>& dims) const;

  /// does this recurrent relation have a generic equivalent? Default is no.
  virtual bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const;
  /// return the name of a header file with the declaration of the generic code
  virtual std::string generic_header() const;
  /// return the implementation of this recurrence relation in terms of generic
  /// code
  virtual std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const;
};

namespace algebra {
/// these operators are extremely useful to write compact expressions
std::shared_ptr<RecurrenceRelation::ExprType> operator+(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B);
std::shared_ptr<RecurrenceRelation::ExprType> operator-(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B);
std::shared_ptr<RecurrenceRelation::ExprType> operator*(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B);
std::shared_ptr<RecurrenceRelation::ExprType> operator/(
    const std::shared_ptr<DGVertex>& A, const std::shared_ptr<DGVertex>& B);
const std::shared_ptr<RecurrenceRelation::ExprType>& operator+=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B);
const std::shared_ptr<RecurrenceRelation::ExprType>& operator-=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B);
const std::shared_ptr<RecurrenceRelation::ExprType>& operator*=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B);
const std::shared_ptr<RecurrenceRelation::ExprType>& operator/=(
    std::shared_ptr<RecurrenceRelation::ExprType>& A,
    const std::shared_ptr<DGVertex>& B);

class Entity;
template <class T>
class RTimeEntity;
template <class T>
class CTimeEntity;
// seems to confound Intel compiler on Linux?
// std::shared_ptr<RecurrenceRelation::ExprType> operator*(const
// std::shared_ptr<Entity>& A,
//                                                const
//                                                std::shared_ptr<DGVertex>& B);
template <typename T>
std::shared_ptr<RecurrenceRelation::ExprType> operator*(
    const std::shared_ptr<RTimeEntity<T> >& A,
    const std::shared_ptr<DGVertex>& B);
template <typename T>
std::shared_ptr<RecurrenceRelation::ExprType> operator*(
    const std::shared_ptr<CTimeEntity<T> >& A,
    const std::shared_ptr<DGVertex>& B);
};  // namespace algebra

// Instantiate the RRStack
typedef RRStackBase<RecurrenceRelation> RRStack;

};  // namespace libint2

// #include <vrr_11_twoprep_11.h>
// #include <hrr.h>
// #include <shell_to_ints.h>
// #include <iter.h>

#endif
