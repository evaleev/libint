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

#ifndef _libint2_src_bin_libint_uncontract_h_
#define _libint2_src_bin_libint_uncontract_h_

#include <algebra.h>
#include <bfset.h>
#include <context.h>
#include <dims.h>
#include <integral.h>
#include <iter.h>
#include <rr.h>

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** Uncontract_Integral_base is dummy class used for dynamic casts only
 */
class Uncontract_Integral_base {
 protected:
  virtual ~Uncontract_Integral_base() {}
};

/** Uncontract_Integral converts (a set of) contracted integral(s) to its
uncontracted counterpart. Although this is technically not a recurrence
relation, it can be treated as one.
*/
template <class I>
class Uncontract_Integral : public RecurrenceRelation,
                            public Uncontract_Integral_base {
 public:
  typedef I TargetType;
  typedef TargetType ChildType;
  typedef typename I::BasisFunctionType BFSet;
  /// The type of expressions in which RecurrenceRelations result.
  typedef RecurrenceRelation::ExprType ExprType;

  Uncontract_Integral(const std::shared_ptr<I>&);
  virtual ~Uncontract_Integral() {}

  /// Implementation of RecurrenceRelation::num_children()
  unsigned int num_children() const override { return children_.size(); }
  /// target() returns pointer to target
  std::shared_ptr<TargetType> target() const { return target_; };
  /// child(i) returns pointer i-th child
  std::shared_ptr<ChildType> child(unsigned int i) const;
  /// Implementation of RecurrenceRelation's target()
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target());
  }
  /// Implementation of RecurrenceRelation's child()
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return std::static_pointer_cast<DGVertex, ChildType>(child(i));
  }
  /// to inline this would require a unary operator (+=).
  /// instead will always implement as a function call.
  bool is_simple() const override { return false; }

 private:
  std::shared_ptr<TargetType> target_;
  std::vector<std::shared_ptr<ChildType> > children_;

  // contracting integrals only depends on the number of integrals in a set
  // can simplify this to only refer to that number
  std::string generate_label() const override {
    std::ostringstream os;
    os << "Generic Contract";
    return os.str();
  }

  //
  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override;
};

template <class I>
Uncontract_Integral<I>::Uncontract_Integral(const std::shared_ptr<I>& Tint)
    : target_(Tint) {
  target_ = Tint;

  // create the uncontracted version
  typedef typename I::BraType bratype;
  typedef typename I::KetType kettype;
  typedef typename I::OperType opertype;
  bratype bra_unc = target_->bra();
  kettype ket_unc = target_->ket();
  // is it even contracted?
  bool target_is_contracted = false;
  {
    const unsigned int np = bra_unc.num_part();
    for (unsigned int p = 0; p < np; ++p) {
      const unsigned int nf = bra_unc.num_members(p);
      for (unsigned int f = 0; f < nf; ++f) {
        target_is_contracted |= bra_unc.member(p, f).contracted();
        bra_unc.member(p, f).uncontract();
      }
    }
  }
  {
    const unsigned int np = ket_unc.num_part();
    for (unsigned int p = 0; p < np; ++p) {
      const unsigned int nf = ket_unc.num_members(p);
      for (unsigned int f = 0; f < nf; ++f) {
        target_is_contracted |= ket_unc.member(p, f).contracted();
        ket_unc.member(p, f).uncontract();
      }
    }
  }
  opertype oper_unc = target_->oper();
  const bool oper_contracted = oper_unc.descr().contracted();
  target_is_contracted |= oper_contracted;
  oper_unc.descr().uncontract();

  if (target_is_contracted) {
#if DEBUG
    std::cout << "Uncontract_Integral: " << target_->description()
              << " is contracted" << std::endl;
#endif
    std::shared_ptr<ChildType> c =
        ChildType::Instance(bra_unc, ket_unc, target_->aux(), oper_unc);
    children_.push_back(c);
  }
};

template <class I>
std::shared_ptr<typename Uncontract_Integral<I>::ChildType>
Uncontract_Integral<I>::child(unsigned int i) const {
  return children_.at(i);
};

template <class I>
std::string Uncontract_Integral<I>::spfunction_call(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims) const {
  const unsigned int s = target_->size();
  std::shared_ptr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));
  std::shared_ptr<Entity> bvecdim;
  bool vectorize = false;
  if (!dims->vecdim_is_static()) {
    vectorize = true;
    std::shared_ptr<RTimeEntity<EntityTypes::Int> > vecdim =
        std::dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>, Entity>(
            dims->vecdim());
    bvecdim = vecdim * bdim;
  } else {
    std::shared_ptr<CTimeEntity<int> > vecdim =
        std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(dims->vecdim());
    vectorize = vecdim->value() == 1 ? false : true;
    bvecdim = vecdim * bdim;
  }

  std::ostringstream os;
  // contraction = reduction
  if (!vectorize || !TrivialBFSet<BFSet>::result ||
      context->cparams()->vectorize_by_line()) {
    os << "_libint2_static_api_inc1_short_("
       << context->value_to_pointer(rr_target()->symbol()) << ","
       << context->value_to_pointer(rr_child(0)->symbol()) << ","
       << bvecdim->id() << ")" << context->end_of_stat() << std::endl;
  } else {  // blockwise vectorize for a single integral
    os << "_libint2_static_api_inc1_short_("
       << context->value_to_pointer(rr_target()->symbol()) << "+vi,"
       << context->value_to_pointer(rr_child(0)->symbol()) << ",1)"
       << context->end_of_stat() << std::endl;
  }
  unsigned int& nflops_ref = const_cast<unsigned int&>(nflops_);
  nflops_ref += target_->size();

  return os.str();
}

/// return true if V is a decontracted IntegralSet
struct DecontractedIntegralSet {
  bool operator()(const std::shared_ptr<DGVertex>& V) {
    const unsigned int outdegree = V->num_exit_arcs();
    if (outdegree == 0) return false;

    const std::shared_ptr<DGArc> arc0 = *(V->first_exit_arc());
    // Is this DGArcRR?
    const std::shared_ptr<DGArcRR> arcrr =
        std::dynamic_pointer_cast<DGArcRR, DGArc>(arc0);
    if (arcrr == 0) return false;
    const std::shared_ptr<Uncontract_Integral_base> uib_ptr =
        std::dynamic_pointer_cast<Uncontract_Integral_base, RecurrenceRelation>(
            arcrr->rr());
    return uib_ptr != 0;
  }
};

};  // namespace libint2

#endif
