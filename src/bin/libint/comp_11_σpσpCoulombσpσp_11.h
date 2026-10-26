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

#ifndef LIBINT_COMP_11_ΣPΣPCOULOMBΣPΣP_11_H
#define LIBINT_COMP_11_ΣPΣPCOULOMBΣPΣP_11_H

#include <dims.h>
#include <entity.h>
#include <gaussoper.h>
#include <generic_rr.h>
#include <rkb_fold_codegen.h>
#include <task.h>
#include <twoprep_11_11.h>

#include <stdexcept>

namespace libint2 {

/**
 * Computes the integral
 * \f$ ((\sigma_1\cdot\hat{p}) a)^\dagger ((\sigma_1\cdot\hat{p}) b)
 *     \frac{1}{r_{12}}
 *     ((\sigma_2\cdot\hat{p}) c)^\dagger ((\sigma_2\cdot\hat{p}) d) \f$
 * over CGShell/CGF by rewriting it as a linear combination of integrals over
 * first derivatives of all four functions (with respect to their centers).
 *
 * The two σ operators act on independent spin spaces (electron 1 on the bra
 * pair a, b; electron 2 on the ket pair c, d). Using the Dirac identity (see
 * e.g. Eq. 1.27 of I. P. Grant, "Relativistic Quantum Theory of Atoms and
 * Molecules", Springer, 2007):
 *   \f$ (\sigma \cdot u)(\sigma \cdot v) = (u \cdot v)I
 *       + i\sigma \cdot (u \times v) \f$
 * applied independently to each electron's spin space gives the tensor
 * product of two 4-component (scalar + vector) objects, with
 * \f$ 4 \times 4 = 16 \f$ components:
 *
 *   index = 4 * q₁ + q₂
 *
 * where q₁ (bra pair) and q₂ (ket pair) are: 0=S (scalar, ∇·∇), 1=X, 2=Y,
 * 3=Z (the components of the cross product ∇×∇). The components are built by
 * fold_terms() as the product of the Pauli components of the two pairs (see
 * rkb_fold_pair_terms):
 *   SS (index 0):                  (∇_a·∇_b)(∇_c·∇_d)            9 terms
 *   SX,SY,SZ (1-3), XS,YS,ZS (4,8,12): dot × cross               6 terms
 *   XX..ZZ (5-7, 9-11, 13-15):     -(∇_a×∇_b)_i (∇_c×∇_d)_j      4 terms
 *
 * Sign convention: the cross × cross block includes the minus sign from
 * \f$ i^2 = -1 \f$, the product of the two \f$ i \f$ factors of the Dirac
 * identity:
 *   \f$ [i\sigma_1 \cdot (u \times v)] \otimes [i\sigma_2 \cdot (w \times x)]
 *     = -\sigma_{1,i} \otimes \sigma_{2,j}\; (u \times v)_i\, (w \times x)_j
 *   \f$
 * so the user multiplies component 4q₁+q₂ by 1 (q₁=q₂=0), i σ_{2,q₂} (q₁=0),
 * i σ_{1,q₁} (q₂=0), or σ_{1,q₁} σ_{2,q₂} (both ≥ 1); the center-derivative
 * convention contributes no further prefactor since
 * (p_i a)† (p_j b) = ∂_{A,i} a ∂_{B,j} b (see the rkb_operators page in
 * engine.h).
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_σpσpCoulombσpσp_11
    : public GenericRecurrenceRelation<
          CR_11_σpσpCoulombσpσp_11<F>, F,
          GenIntegralSet_11_11<F, σpσpCoulombσpσpOper, mType>> {
 public:
  typedef CR_11_σpσpCoulombσpσp_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpσpCoulombσpσpOper OperType;
  typedef GenIntegralSet_11_11<F, σpσpCoulombσpσpOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  // children_.reserve() hint only (see generic_rr.h): the SS component sums
  // 9 deriv-ERI children, all others 6 or 4
  static const unsigned int max_nchildren = 9;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by Instance that maintains
  /// registry of these objects
  CR_11_σpσpCoulombσpσp_11(const std::shared_ptr<TargetType>&,
                           unsigned int = 0);

  static std::string descr() { return "CR"; }

  /// The terms of component @p comp = 4 q₁ + q₂: the product of the terms of
  /// Pauli component q₁ of (σ₁·∇_a)(σ₁·∇_b) (derivative directions on a = bra1
  /// and b = ket1) and of q₂ of (σ₂·∇_c)(σ₂·∇_d) (on c = bra2 and d = ket2),
  /// see rkb_fold_pair_terms; the i's of two cross products multiply to -1.
  /// The constructor (the children) and generate_code() (the element-wise
  /// expression) are both driven by this table, so they cannot disagree.
  static RKBFoldTerms fold_terms(int comp) {
    if (comp < 0 || comp > 15)
      throw std::runtime_error(
          "CR_11_σpσpCoulombσpσp_11: invalid component index " +
          std::to_string(comp) + " (expected 0-15)");
    const int q1 = comp / 4;
    const int q2 = comp % 4;
    const int sign = (q1 != 0 && q2 != 0) ? -1 : 1;
    RKBFoldTerms terms;
    for (const auto& p1 : rkb_fold_pair_terms(q1))
      for (const auto& p2 : rkb_fold_pair_terms(q2))
        terms.push_back(
            RKBFoldTerm{sign * p1.coef * p2.coef, {{p1.i, p1.j, p2.i, p2.j}}});
    return terms;
  }

  // --- Code sharing overrides ---
  // All shell sets with the same component share one function.

  std::string generate_label() const override {
    return "CR_opopCoulombopop_" +
           std::to_string(target_->oper()->descr().component_index());
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
  /// target = src0 + src1 + ... + src8 for the SS component
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    rkb_fold_generate_code_from_terms(
        *this, fold_terms(target_->oper()->descr().component_index()),
        "CR_11_σpσpCoulombσpσp_11", context, dims, funcname, decl, def);
  }
};

template <typename F>
CR_11_σpσpCoulombσpσp_11<F>::CR_11_σpσpCoulombσpσp_11(
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
                         fold_terms(Tint->oper()->descr().component_index()),
                         is_simple(), expr_, nflops_);

}  // CR_11_σpσpCoulombσpσp_11<F>::CR_11_σpσpCoulombσpσp_11

}  // namespace libint2

#endif  // LIBINT_COMP_11_ΣPΣPCOULOMBΣPΣP_11_H
