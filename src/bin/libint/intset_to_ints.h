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
#include <integral.h>
#include <iter.h>
#include <rr.h>

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef _libint2_src_bin_libint_intsettoints_h_
#define _libint2_src_bin_libint_intsettoints_h_

namespace libint2 {

/** IntegralSet_to_Integrals_base is dummy class used for dynamic casts only
 */
class IntegralSet_to_Integrals_base {
 protected:
  virtual ~IntegralSet_to_Integrals_base() {}
};

/** IntegralSet_to_Integrals converts I, a set of integrals, to individual
integrals. Although this is technically not a recurrence relation, it can be
treated as one.
*/
template <class I>
class IntegralSet_to_Integrals : public RecurrenceRelation,
                                 public IntegralSet_to_Integrals_base {
 public:
  typedef I TargetType;
  typedef typename I::iter_type ChildType;
  /// The type of expressions in which RecurrenceRelations result.
  typedef RecurrenceRelation::ExprType ExprType;

  IntegralSet_to_Integrals(const std::shared_ptr<I>&);
  virtual ~IntegralSet_to_Integrals() {}

  /// Return an instance if applicable, or a null pointer otherwise
  static std::shared_ptr<IntegralSet_to_Integrals<I>> Instance(
      const std::shared_ptr<TargetType>& Tint, unsigned int dir) {
    assert(dir == 0);
    // attempt to construct
    std::shared_ptr<IntegralSet_to_Integrals<I>> this_ptr(
        new IntegralSet_to_Integrals<I>(Tint));
    // if succeeded (nchildren > 0) do post-construction
    assert(this_ptr->num_children() != 0);
    return this_ptr;
  }
  static bool directional() { return false; }

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
  /// Implementation of RecurrenceRelation::is_simple()
  bool is_simple() const override { return true; }
  /// Reimplementation of RecurrenceRelation::invariant_type()
  bool invariant_type() const override {
    // Converts from one BFSet to another!
    return false;
  }

 private:
  std::shared_ptr<TargetType> target_;
  std::vector<std::shared_ptr<ChildType>> children_;

  /// Implementation of RecurrenceRelation::generate_label()
  std::string generate_label() const override {
    return "IntegralSet_to_Integrals";
    // throw std::runtime_error("IntegralSet_to_Integrals::label() -- code for
    // this RR is never generated, so this function should never be used");
  }
  /// Reimplementation of RecurrenceRelation::spfunction_call()
  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override {
    throw std::logic_error(
        "IntegralSet_to_Integrals::spfunction_call -- should not call this "
        "function");
  }
};

template <class I>
IntegralSet_to_Integrals<I>::IntegralSet_to_Integrals(
    const std::shared_ptr<I>& Tint)
    : target_(Tint) {
  target_ = Tint;

  // Construct a subiterator for I
  SubIteratorBase<I> siter(Tint);

  // Set children pointers
  for (siter.init(); siter; ++siter) children_.push_back(siter.elem());
};

template <class I>
std::shared_ptr<typename I::iter_type> IntegralSet_to_Integrals<I>::child(
    unsigned int i) const {
  return children_.at(i);
};

};  // namespace libint2

#endif
