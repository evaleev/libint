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

#ifndef _libint2_src_bin_libint_cr11r2dotr2g1211_h_
#define _libint2_src_bin_libint_cr11r2dotr2g1211_h_

#include <algebra.h>
#include <context.h>
#include <default_params.h>
#include <dgvertex.h>
#include <flop.h>
#include <integral.h>
#include <prefactors.h>
#include <r2dotr2g12_11_11.h>
#include <rr.h>
#include <rr.templ.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** Compute relation for 2-e integrals of the r2.r2 x G12 operators.
I<BFSet> is the integral set specialization that describes the
integrals of the R2dotR2_G12 operator.
*/
template <template <class> class I, class BFSet>
class CR_11_R2dotR2G12_11 : public RecurrenceRelation {
 public:
  typedef RecurrenceRelation ParentType;
  typedef BFSet BasisFunctionType;
  typedef CR_11_R2dotR2G12_11<I, BFSet> ThisType;
  typedef I<BFSet> TargetType;
  typedef R12kG12_11_11<BFSet, 0> ChildType;
  /// The type of expressions in which RecurrenceRelations result.
  typedef RecurrenceRelation::ExprType ExprType;

  /** Use Instance() to obtain an instance of RR. This function is provided to
     avoid issues with getting a std::shared_ptr from constructor (as needed for
     registry to work).
  */
  static std::shared_ptr<ThisType> Instance(const std::shared_ptr<TargetType>&);
  virtual ~CR_11_R2dotR2G12_11() {}

  /// Implementation of RecurrenceRelation::num_children()
  unsigned int num_children() const override { return nchildren_; }
  /// target() returns pointer to the i-th child
  std::shared_ptr<TargetType> target() const { return target_; };
  /// child(i) returns pointer to the i-th child
  std::shared_ptr<ChildType> child(unsigned int i) const;
  /// Implementation of RecurrenceRelation::rr_target()
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target());
  }
  /// Implementation of RecurrenceRelation::rr_child()
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return std::dynamic_pointer_cast<DGVertex, ChildType>(child(i));
  }
  /// Implementation of RecurrenceRelation::is_simple()
  bool is_simple() const override { return TrivialBFSet<BFSet>::result; }

 private:
  /**
    dir specifies which quantum number is incremented.
    For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
    a Cartesian Gaussian.
   */
  CR_11_R2dotR2G12_11(const std::shared_ptr<TargetType>&);

#if 0
    /// registers this RR with the stack, if needed
    bool register_with_rrstack() const;
#endif

  static const unsigned int max_nchildren_ = 3;
  std::shared_ptr<TargetType> target_;
  std::shared_ptr<ChildType> children_[max_nchildren_];
  unsigned int nchildren_;

  std::string generate_label() const override {
    std::ostringstream os;
    os << "RR ( " << rr_target()->label() << " )";
    return os.str();
  }
};

template <template <class> class I, class F>
std::shared_ptr<CR_11_R2dotR2G12_11<I, F> > CR_11_R2dotR2G12_11<I, F>::Instance(
    const std::shared_ptr<TargetType>& Tint) {
  std::shared_ptr<ThisType> this_ptr(new ThisType(Tint));
  // Do post-construction duties
  if (this_ptr->num_children() != 0) {
    this_ptr->register_with_rrstack<ThisType>();
    return this_ptr;
  }
  return std::shared_ptr<ThisType>();
}

template <template <class> class I, class F>
CR_11_R2dotR2G12_11<I, F>::CR_11_R2dotR2G12_11(
    const std::shared_ptr<I<F> >& Tint)
    : ParentType(), target_(Tint), nchildren_(0) {
  F sh_a(Tint->bra(0, 0));
  F sh_b(Tint->ket(0, 0));
  F sh_c(Tint->bra(1, 0));
  F sh_d(Tint->ket(1, 0));

  std::vector<F> bra;
  std::vector<F> ket;
  bra.push_back(sh_a);
  bra.push_back(sh_c);
  ket.push_back(sh_b);
  ket.push_back(sh_d);

  auto* bra_ref = &bra;
  auto* ket_ref = &ket;

  const unsigned int ndirs = is_simple() ? 3 : 1;
  for (int xyz = 0; xyz < ndirs; xyz++) {
    // c+1_i d+1_i
    bra_ref->operator[](1).inc(xyz);
    ket_ref->operator[](1).inc(xyz);
    int next_child = nchildren_;
    children_[next_child] =
        ChildType::Instance(bra[0], ket[0], bra[1], ket[1], 0);
    ++nchildren_;

    if (is_simple()) {
      std::shared_ptr<ExprType> expr0_ptr(new ExprType(
          ExprType::OperatorTypes::Times, Scalar(1.0), rr_child(next_child)));
      add_expr(expr0_ptr);
      nflops_ += 1;
    }
    bra_ref->operator[](1).dec(xyz);
    ket_ref->operator[](1).dec(xyz);
  }
}

template <template <class> class I, class F>
std::shared_ptr<typename CR_11_R2dotR2G12_11<I, F>::ChildType>
CR_11_R2dotR2G12_11<I, F>::child(unsigned int i) const {
  assert(i >= 0 && i < nchildren_);
  unsigned int nc = 0;
  for (int c = 0; c < max_nchildren_; c++) {
    if (children_[c] != 0) {
      if (nc == i) return children_[c];
      nc++;
    }
  }
};

/// Useful typedefs
typedef CR_11_R2dotR2G12_11<R2dotR2G12_11_11, CGShell> CR_11_R2dotR2G12_11_sq;
typedef CR_11_R2dotR2G12_11<R2dotR2G12_11_11, CGF> CR_11_R2dotR2G12_11_int;

};  // namespace libint2

#endif
