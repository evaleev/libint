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

#ifndef LIBINT_COMP_11_COULOMBΣPΣP_11_H
#define LIBINT_COMP_11_COULOMBΣPΣP_11_H

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
 * Computes the integral (a b | 1/r_{12} | (σ·p c)† (σ·p d)) over CGShell/CGF,
 * with σ·p acting on both functions of electron 2 (c = bra2, d = ket2), as a
 * linear combination of first-derivative TwoPRep integrals. The Dirac identity
 * σ_i σ_j = δ_ij + i ε_ijk σ_k folds the two Pauli matrices into 4 quaternion
 * components: 0: ∇_c·∇_d, 1-3: (∇_c × ∇_d)_{x,y,z}. The derivatives are with
 * respect to the centers of c and d; since (p_i c)† (p_j d) =
 * (-i ∂_{C,i} c)(+i ∂_{D,j} d) = ∂_{C,i} c ∂_{D,j} d, the user assembles the
 * 2×2 spin block as q₀ + i (q₁ σ_x + q₂ σ_y + q₃ σ_z) with no further
 * prefactor (see the rkb_operators page in engine.h).
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_Coulombσpσp_11
    : public GenericRecurrenceRelation<
          CR_11_Coulombσpσp_11<F>, F,
          GenIntegralSet_11_11<F, CoulombσpσpOper, mType>> {
 public:
  typedef CR_11_Coulombσpσp_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef CoulombσpσpOper OperType;
  typedef GenIntegralSet_11_11<F, CoulombσpσpOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  // children_.reserve() hint only (see generic_rr.h): the scalar component
  // sums 3 deriv-ERI children, the cross-product components 2
  static const unsigned int max_nchildren = 3;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by Instance that maintains
  /// registry of these objects
  CR_11_Coulombσpσp_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  /// The terms of quaternion component @p q of (σ·∇_c)(σ·∇_d) (see
  /// rkb_fold_pair_terms): derivative direction i on c (bra2), j on d (ket2).
  /// The constructor (the children) and generate_code() (the element-wise
  /// expression) are both driven by this table, so they cannot disagree.
  static RKBFoldTerms fold_terms(int q) {
    RKBFoldTerms terms;
    for (const auto& p : rkb_fold_pair_terms(q))
      terms.push_back(RKBFoldTerm{p.coef, {{-1, -1, p.i, p.j}}});
    return terms;
  }

  // --- Code sharing overrides ---
  // All shell sets with the same quaternion component share one function.

  std::string generate_label() const override {
    return "CR_Coulombopop_" +
           std::to_string(target_->oper()->descr().quaternion_index());
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

  /// emits the element-wise fold of the children, e.g.
  /// target = src0 + src1 + src2 for the scalar component
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    rkb_fold_generate_code_from_terms(
        *this, fold_terms(target_->oper()->descr().quaternion_index()),
        "CR_11_Coulombσpσp_11", context, dims, funcname, decl, def);
  }
};

template <typename F>
CR_11_Coulombσpσp_11<F>::CR_11_Coulombσpσp_11(
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
                         fold_terms(Tint->oper()->descr().quaternion_index()),
                         is_simple(), expr_, nflops_);

}  // CR_11_Coulombσpσp_11<F>::CR_11_Coulombσpσp_11

}  // namespace libint2

#endif  // LIBINT_COMP_11_COULOMBΣPΣP_11_H
