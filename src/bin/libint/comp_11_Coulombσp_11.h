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

#ifndef LIBINT_COMP_11_COULOMBΣP_11_H
#define LIBINT_COMP_11_COULOMBΣP_11_H

#include <dims.h>
#include <entity.h>
#include <gaussoper.h>
#include <generic_rr.h>
#include <task.h>
#include <twoprep_11_11.h>

#include <rkb_fold_codegen.h>

namespace libint2 {

/**
 * this computes the 3-center single-σ·p integral
 * \f$ \frac{1}{r_{12}} \, \sigma \cdot \hat{p}_2 \f$ acting on the 2nd ket
 * function ν, by rewriting it as the bare Cartesian gradient of ν, i.e. one
 * first-derivative integral over \f$ \frac{1}{r_{12}} \f$ per Cartesian
 * direction. Unlike Coulombσpσp there is no Dirac σ·σ folding: each of the 3
 * components (x,y,z) is exactly one derivative ERI child, so it routes
 * naturally through the TwoPRep DerivGaussV2 strategy.
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_Coulombσp_11
    : public GenericRecurrenceRelation<
          CR_11_Coulombσp_11<F>, F,
          GenIntegralSet_11_11<F, CoulombσpOper, mType>> {
 public:
  typedef CR_11_Coulombσp_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef CoulombσpOper OperType;
  typedef GenIntegralSet_11_11<F, CoulombσpOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 1;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by Instance that maintains
  /// registry of these objects
  CR_11_Coulombσp_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  // --- Code sharing overrides ---
  // All shell combos with the same Cartesian component share one function
  // (the element-wise copy below is spectator-dimension-independent).

  std::string generate_label() const override {
    return "CR_Coulombop_" +
           std::to_string(target_->oper()->descr().cartesian_index());
  }

  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override {
    return rkb_fold_spfunction_call(*this, target_, context);
  }

  std::shared_ptr<ImplicitDimensions> adapt_dims_(
      const std::shared_ptr<ImplicitDimensions>& dims) const override {
    return rkb_fold_adapt_dims(dims);
  }

  /// Hand-emit a simple element-wise copy: target = src0 (the single
  /// first-derivative child). 0 FLOPs per component.
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    // target = src0 (the single first-derivative child), 0 flops.
    rkb_fold_generate_code(*this, context, dims, funcname, decl, def,
                           "src0[hsi]", 0);
  }
};

template <typename F>
CR_11_Coulombσp_11<F>::CR_11_Coulombσp_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_bra(/* particle */ 1) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 1) == 1);

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  const auto& oper = Tint->oper();

  if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
    return;

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_11_11<BasisFunctionType, TwoPRep, mType>>
      factory(this);

  // single σ·p on ν (the 2nd ket function, d). Component k = ∂(d)/∂{x,y,z}.
  // (Geometric-center derivative convention, matching CR_11_Coulombσpσp_11; the
  // -i factor and Pauli matrix are applied by the caller.)
  const int k = oper->descr().cartesian_index();
  if (k < 0 || k > 2)
    throw std::runtime_error("CR_11_Coulombσp_11: invalid Cartesian index");

  F d_k{d};
  d_k.deriv().inc(k);  // d(d)/d{x,y,z} = d_k

  auto a_b_c_dk = factory.make_child(a, b, c, d_k, zero_m);
  if (is_simple()) {
    // Single derivative child; Scalar(1)* is only a registration vehicle (a
    // bare child is not an expression). The actual emit is the element-wise
    // copy in generate_code() above, so this contributes 0 real flops.
    expr_ = Scalar(1) * a_b_c_dk;
    nflops_ += 0;
  }

}  // CR_11_Coulombσp_11<F>::CR_11_Coulombσp_11
};  // namespace libint2

#endif  // LIBINT_COMP_11_COULOMBΣP_11_H
