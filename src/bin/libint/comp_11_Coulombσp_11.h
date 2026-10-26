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

#include <stdexcept>

namespace libint2 {

/**
 * Computes the 3-center single-σ·p integral (a b | 1/r_{12} | c σ·p d), with
 * the fitting function as a (bra1) and the unit function as b (ket1). Each
 * Cartesian component k is the TwoPRep integral over d differentiated with
 * respect to its center D in direction k (libint's derivative convention),
 * i.e. exactly one first-derivative ERI child, so it routes naturally through
 * the TwoPRep DerivGaussV2 strategy; unlike Coulombσpσp there is no σ·σ
 * folding. Since p_k = -i ∂/∂r_k = +i ∂/∂D_k, the user multiplies component k
 * by +i and contracts with σ_k (see the rkb_operators page in engine.h).
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_Coulombσp_11 : public GenericRecurrenceRelation<
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

  /// The single term of Cartesian component @p k: the first derivative of d
  /// (ket2) with respect to its center, direction k. The constructor (the
  /// children) and generate_code() (the element-wise expression) are both
  /// driven by this table, so they cannot disagree.
  static RKBFoldTerms fold_terms(int k) {
    if (k < 0 || k > 2)
      throw std::runtime_error("CR_11_Coulombσp_11: invalid Cartesian index " +
                               std::to_string(k));
    return {RKBFoldTerm{1, {{-1, -1, -1, k}}}};
  }

  // --- Code sharing overrides ---
  // All shell sets with the same Cartesian component share one function
  // (the element-wise copy is spectator-dimension-independent).

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

  /// emits the element-wise copy target = src0 (0 flops)
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    rkb_fold_generate_code_from_terms(
        *this, fold_terms(target_->oper()->descr().cartesian_index()),
        "CR_11_Coulombσp_11", context, dims, funcname, decl, def);
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

  if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
    return;

  ChildFactory<ThisType,
               GenIntegralSet_11_11<BasisFunctionType, TwoPRep, mType>>
      factory(this);
  rkb_fold_make_children(factory, a, b, c, d,
                         fold_terms(Tint->oper()->descr().cartesian_index()),
                         is_simple(), expr_, nflops_);

}  // CR_11_Coulombσp_11<F>::CR_11_Coulombσp_11

}  // namespace libint2

#endif  // LIBINT_COMP_11_COULOMBΣP_11_H
