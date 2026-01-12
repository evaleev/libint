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

#ifndef _libint2_src_bin_libint_itr11twoprep11_h_
#define _libint2_src_bin_libint_itr11twoprep11_h_

#include <algebra.h>
#include <context.h>
#include <default_params.h>
#include <dgvertex.h>
#include <flop.h>
#include <prefactors.h>
#include <rr.h>
#include <twoprep_11_11.h>
#include <util.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** ITR (Interelectron Transfer Relation) for 2-e ERI. part specifies for which
particle the angular momentum is raised. where specifies whether the angular
momentum is shifted in bra or ket. Class ERI specifies which particular
implementation of ERI to use.
*/
template <template <typename, typename, typename> class ERI, class BFSet,
          int part, FunctionPosition where>
class ITR_11_TwoPRep_11 : public RecurrenceRelation {
 public:
  typedef RecurrenceRelation ParentType;
  typedef BFSet BasisFunctionType;
  typedef ITR_11_TwoPRep_11 ThisType;
  typedef ERI<BFSet, TwoPRep, mType> TargetType;
  typedef TargetType ChildType;
  /// The type of expressions in which RecurrenceRelations result.
  typedef RecurrenceRelation::ExprType ExprType;

  /** Use Instance() to obtain an instance of RR. This function is provided to
     avoid issues with getting a std::shared_ptr from constructor (as needed for
     registry to work).

      dir specifies which quantum number of a and b is shifted.
      For example, dir can be 0 (x), 1(y), or 2(z) if F is
      a Cartesian Gaussian.
  */
  static std::shared_ptr<ThisType> Instance(const std::shared_ptr<TargetType>&,
                                            unsigned int dir = 0);
  virtual ~ITR_11_TwoPRep_11() { assert(part == 0 || part == 1); }

  /// overrides RecurrenceRelation::partindex_direction()
  int partindex_direction() const {
    return part == 0 ? +1   // transfer from 0 to 1
                     : -1;  // transfer from 1 to 0
  }

  /// Default directionality
  /** is this recurrence relation parameterized by a direction.
      the default is false if BasisFunctionSet is CGShell,
      true otherwise. */
  static bool directional() {
    if (boost::is_same<BasisFunctionType, CGShell>::value) return false;
    return true;
  }

#if 1
  /// Implementation of RecurrenceRelation::num_children()
  unsigned int num_children() const override { return children_.size(); }
  /// Implementation of RecurrenceRelation::rr_target()
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target_);
  }
  /// Implementation of RecurrenceRelation::rr_child()
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return std::static_pointer_cast<DGVertex, ChildType>(children_.at(i));
  }
  /// Implementation of RecurrenceRelation::is_simple()
  bool is_simple() const override { return TrivialBFSet<BFSet>::result; }
#endif

 private:
  /**
    dir specifies which quantum number is incremented.
    For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
    a Cartesian Gaussian.
   */
  ITR_11_TwoPRep_11(const std::shared_ptr<TargetType>&, unsigned int dir);

  static const unsigned int max_nchildren_ = 4;
  unsigned int dir_;
  std::shared_ptr<TargetType> target_;
  std::vector<std::shared_ptr<ChildType> > children_;
  const std::shared_ptr<ChildType>& make_child(const BFSet& A, const BFSet& B,
                                               const BFSet& C, const BFSet& D,
                                               unsigned int m) {
    const std::shared_ptr<ChildType>& i = ChildType::Instance(A, B, C, D, m);
    children_.push_back(i);
    return *(children_.end() - 1);
  }

  std::string generate_label() const override {
    typedef typename TargetType::AuxIndexType mType;
    static std::shared_ptr<mType> aux0(new mType(0u));
    std::ostringstream os;
    // ITR recurrence relation code is independent of m (it never appears
    // anywhere in equations), hence to avoid generating identical code make
    // sure that the (unique) label does not contain m
    os << "ITR Part" << part << " " << to_string(where)
       << genintegralset_label(target_->bra(), target_->ket(), aux0,
                               target_->oper());
    return os.str();
  }

#if LIBINT_ENABLE_GENERIC_CODE
  /// Implementation of RecurrenceRelation::has_generic()
  bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const override;
  /// Implementation of RecurrenceRelation::generic_header()
  std::string generic_header() const override;
  /// Implementation of RecurrenceRelation::generic_instance()
  std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const override;
#endif
};

template <template <typename, typename, typename> class ERI, class F, int part,
          FunctionPosition where>
std::shared_ptr<ITR_11_TwoPRep_11<ERI, F, part, where> >
ITR_11_TwoPRep_11<ERI, F, part, where>::Instance(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir) {
  std::shared_ptr<ThisType> this_ptr(new ThisType(Tint, dir));
  // Do post-construction duties
  if (this_ptr->num_children() != 0) {
    this_ptr->template register_with_rrstack<ThisType>();
    return this_ptr;
  }
  return std::shared_ptr<ThisType>();
}

template <template <typename, typename, typename> class ERI, class F, int part,
          FunctionPosition where>
ITR_11_TwoPRep_11<ERI, F, part, where>::ITR_11_TwoPRep_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : target_(Tint), dir_(dir) {
  // only works for primitive integrals
  {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));
    if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
      return;
  }

  // not yet implemented for derivative integrals
  {
    const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
    const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
    const OriginDerivative<3u> dC = Tint->bra(1, 0).deriv();
    const OriginDerivative<3u> dD = Tint->ket(1, 0).deriv();
    const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();
    if (deriv) return;
  }

  children_.reserve(max_nchildren_);
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  const unsigned int m = Tint->aux()->elem(0);
  const F& _1 = unit<F>(dir);

  // Build on A
  if (part == 0 && where == InBra) {
    F a(Tint->bra(0, 0) - _1);
    if (!exists(a)) return;
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));

    // must be (A0|C0)
    if (not(b.zero() && d.zero())) return;

    const std::shared_ptr<ChildType>& ABCD_m = make_child(a, b, c, d, m);
    if (is_simple()) {
      expr_ = Vector("TwoPRepITR_pfac0_0_0")[dir] * ABCD_m;
      nflops_ += 1;
    }

    const F& cp1 = c + _1;
    const std::shared_ptr<ChildType>& ABCp1D_m = make_child(a, b, cp1, d, m);
    if (is_simple()) {
      expr_ -= Scalar("eoz") * ABCp1D_m;
      nflops_ += 2;
    }

    const F& am1 = a - _1;
    if (exists(am1)) {
      const std::shared_ptr<ChildType>& Am1BCD_m = make_child(am1, b, c, d, m);
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2z") * Am1BCD_m;
        nflops_ += 3;
      }
    }
    const F& cm1 = c - _1;
    if (exists(cm1)) {
      const std::shared_ptr<ChildType>& ABCm1D_m = make_child(a, b, cm1, d, m);
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2z") * ABCm1D_m;
        nflops_ += 3;
      }
    }
    return;
  }
  // Build on C
  if (part == 1 && where == InBra) {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0) - _1);
    if (!exists(c)) return;
    F d(Tint->ket(1, 0));

    // must be (A0|C0)
    if (not(b.zero() && d.zero())) return;

    const std::shared_ptr<ChildType>& ABCD_m = make_child(a, b, c, d, m);
    if (is_simple()) {
      expr_ = Vector("TwoPRepITR_pfac0_1_0")[dir] * ABCD_m;
      nflops_ += 1;
    }

    const F& ap1 = a + _1;
    const std::shared_ptr<ChildType>& Ap1BCD_m = make_child(ap1, b, c, d, m);
    if (is_simple()) {
      expr_ -= Scalar("zoe") * Ap1BCD_m;
      nflops_ += 2;
    }

    const F& cm1 = c - _1;
    if (exists(cm1)) {
      const std::shared_ptr<ChildType>& ABCm1D_m = make_child(a, b, cm1, d, m);
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2e") * ABCm1D_m;
        nflops_ += 3;
      }
    }
    const F& am1 = a - _1;
    if (exists(am1)) {
      const std::shared_ptr<ChildType>& Am1BCD_m = make_child(am1, b, c, d, m);
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2e") * Am1BCD_m;
        nflops_ += 3;
      }
    }
    return;
  }
  // Build on B
  if (part == 0 && where == InKet) {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0) - _1);
    if (!exists(b)) return;
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));

    // must be (0B|0D)
    if (not(a.zero() && c.zero())) return;

    const std::shared_ptr<ChildType>& ABCD_m = make_child(a, b, c, d, m);
    if (is_simple()) {
      expr_ = Vector("TwoPRepITR_pfac0_0_1")[dir] * ABCD_m;
      nflops_ += 1;
    }

    const F& dp1 = d + _1;
    const std::shared_ptr<ChildType>& ABCDp1_m = make_child(a, b, c, dp1, m);
    if (is_simple()) {
      expr_ -= Scalar("eoz") * ABCDp1_m;
      nflops_ += 2;
    }

    const F& bm1 = b - _1;
    if (exists(bm1)) {
      const std::shared_ptr<ChildType>& ABm1CD_m = make_child(a, bm1, c, d, m);
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2z") * ABm1CD_m;
        nflops_ += 3;
      }
    }
    const F& dm1 = d - _1;
    if (exists(dm1)) {
      const std::shared_ptr<ChildType>& ABCDm1_m = make_child(a, b, c, dm1, m);
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2z") * ABCDm1_m;
        nflops_ += 3;
      }
    }
    return;
  }
  // Build on D
  if (part == 1 && where == InKet) {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0) - _1);
    if (!exists(d)) return;

    // must be (0B|0D)
    if (not(a.zero() && c.zero())) return;

    const std::shared_ptr<ChildType>& ABCD_m = make_child(a, b, c, d, m);
    if (is_simple()) {
      expr_ = Vector("TwoPRepITR_pfac0_1_1")[dir] * ABCD_m;
      nflops_ += 1;
    }

    const F& bp1 = b + _1;
    const std::shared_ptr<ChildType>& ABp1CD_m = make_child(a, bp1, c, d, m);
    if (is_simple()) {
      expr_ -= Scalar("zoe") * ABp1CD_m;
      nflops_ += 2;
    }

    const F& dm1 = d - _1;
    if (exists(dm1)) {
      const std::shared_ptr<ChildType>& ABCDm1_m = make_child(a, b, c, dm1, m);
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2e") * ABCDm1_m;
        nflops_ += 3;
      }
    }
    const F& bm1 = b - _1;
    if (exists(bm1)) {
      const std::shared_ptr<ChildType>& ABm1CD_m = make_child(a, bm1, c, d, m);
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2e") * ABm1CD_m;
        nflops_ += 3;
      }
    }
    return;
  }
}

#if LIBINT_ENABLE_GENERIC_CODE
template <template <typename, typename, typename> class ERI, class F, int part,
          FunctionPosition where>
bool ITR_11_TwoPRep_11<ERI, F, part, where>::has_generic(
    const std::shared_ptr<CompilationParameters>& cparams) const {
  if (TrivialBFSet<F>::result) return false;

  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const unsigned int max_opt_am = cparams->max_am_opt();
  // generic code works for a0c0 of 0a0c classes where am(a) > 1 and am(c) > 1
  // to generate optimized code for xxxx integral need to generate specialized
  // code for up to (x+x)0(x+x)0 integrals
  if (sh_b.zero() && sh_d.zero() &&
      (sh_a.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_c.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_a.norm() > 1u && sh_c.norm() > 1u)) {
    const bool ITR_xs_xs_Part1_implemented =
        false;  // only Part0 is implemented
    if (part == 1)
      return ITR_xs_xs_Part1_implemented;
    else
      return true;
  }
  if (sh_a.zero() && sh_c.zero() &&
      (sh_b.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_d.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_b.norm() > 1u && sh_d.norm() > 1u)) {
    const bool ITR_sx_sx_implemented = false;  // only ITR_xs_xs is implemented
    return ITR_sx_sx_implemented;
  }
  return false;
}

template <template <typename, typename, typename> class ERI, class F, int part,
          FunctionPosition where>
std::string ITR_11_TwoPRep_11<ERI, F, part, where>::generic_header() const {
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();

  const OriginDerivative<3u> dA = target_->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = target_->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = target_->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = target_->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();
  assert(!deriv);

  if (!deriv) {
    return std::string("ITR_xs_xs.h");
  } else {
    if (xsxs) return std::string("ITR_xs_xs_deriv.h");
    if (sxsx) return std::string("ITR_sx_sx_deriv.h");
  }
  abort();  // unreachable
}

template <template <typename, typename, typename> class ERI, class F, int part,
          FunctionPosition where>
std::string ITR_11_TwoPRep_11<ERI, F, part, where>::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  std::ostringstream oss;
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();

  const OriginDerivative<3u> dA = target_->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = target_->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = target_->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = target_->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();
  assert(!deriv);

  oss << "using namespace libint2;" << std::endl;

  if (!deriv) {  // for regular integrals I know exactly how many prerequisites
                 // I need
    if (xsxs) {
      oss << "libint2::ITR_xs_xs<" << part << "," << sh_a.norm() << ","
          << sh_c.norm() << ",true,";
    }
    if (sxsx) {
      oss << "libint2::ITR_xs_xs<" << part << "," << sh_b.norm() << ","
          << sh_d.norm() << ",false,";
    }
    oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
    oss << ">::compute(inteval";

    const unsigned int nargs = args->n();
    for (unsigned int a = 0; a < nargs; a++) {
      oss << "," << args->symbol(a);
    }
    oss << ");";
  }

  return oss.str();
}
#endif  // #if !LIBINT_ENABLE_GENERIC_CODE

};  // namespace libint2

#endif
