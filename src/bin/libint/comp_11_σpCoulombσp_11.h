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

#ifndef LIBINT_COMP_11_ΣPCOULOMBΣP_11_H
#define LIBINT_COMP_11_ΣPCOULOMBΣP_11_H

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
 * Computes the "Gaunt LS bilinear" integral (a σ·p b | 1/r_{12} | c σ·p d),
 * with one σ·p on the second function of each electron (b = ket1 in the spin
 * space of electron 1, d = ket2 in that of electron 2), in the SO(3)
 * irreducible decomposition of the rank-2 tensor
 * T_ij = (a ∂_{B,i} b | 1/r_{12} | c ∂_{D,j} d): 1 scalar trace + 3
 * antisymmetric (curl-curl) + 5 symmetric-traceless = 9 components, each a
 * small linear combination of first-derivative TwoPRep integrals. The
 * derivatives are with respect to the centers B and D; since
 * (p_i b)(p_j d) = (+i ∂_{B,i} b)(+i ∂_{D,j} d) = -∂_{B,i} b ∂_{D,j} d, the
 * user multiplies the tensor by -1 before contracting with σ_{1,i} σ_{2,j}
 * (see the rkb_operators page in engine.h).
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_σpCoulombσp_11
    : public GenericRecurrenceRelation<
          CR_11_σpCoulombσp_11<F>, F,
          GenIntegralSet_11_11<F, σpCoulombσpOper, mType>> {
 public:
  typedef CR_11_σpCoulombσp_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpCoulombσpOper OperType;
  typedef GenIntegralSet_11_11<F, σpCoulombσpOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
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
  CR_11_σpCoulombσp_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  /// The terms of SO(3) irreducible component @p comp of T_ij: coefficient,
  /// derivative direction i on b (ket1) and j on d (ket2). The constructor
  /// (the children) and generate_code() (the element-wise expression) are both
  /// driven by this table, so they cannot disagree.
  static RKBFoldTerms fold_terms(int comp) {
    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;
    auto T = [](int coef, int i, int j) {
      return RKBFoldTerm{coef, {{-1, i, -1, j}}};
    };
    switch (comp) {
      case σpCoulombσp_Descr::Scalar:  // T_xx + T_yy + T_zz
        return {T(1, x, x), T(1, y, y), T(1, z, z)};
      case σpCoulombσp_Descr::AntisymX:  // T_yz - T_zy
        return {T(1, y, z), T(-1, z, y)};
      case σpCoulombσp_Descr::AntisymY:  // T_zx - T_xz
        return {T(1, z, x), T(-1, x, z)};
      case σpCoulombσp_Descr::AntisymZ:  // T_xy - T_yx
        return {T(1, x, y), T(-1, y, x)};
      case σpCoulombσp_Descr::SymTLDiagA:  // T_xx - T_yy
        return {T(1, x, x), T(-1, y, y)};
      case σpCoulombσp_Descr::SymTLDiagB:  // 2 T_zz - T_xx - T_yy
        return {T(2, z, z), T(-1, x, x), T(-1, y, y)};
      case σpCoulombσp_Descr::SymTLOffXY:  // T_xy + T_yx
        return {T(1, x, y), T(1, y, x)};
      case σpCoulombσp_Descr::SymTLOffXZ:  // T_xz + T_zx
        return {T(1, x, z), T(1, z, x)};
      case σpCoulombσp_Descr::SymTLOffYZ:  // T_yz + T_zy
        return {T(1, y, z), T(1, z, y)};
      default:
        throw std::runtime_error(
            "CR_11_σpCoulombσp_11: invalid component index " +
            std::to_string(comp));
    }
  }

  // --- Code sharing overrides (mirror Coulombσpσp pattern) ---
  // All shell sets with the same component share one function.

  std::string generate_label() const override {
    return "CR_opCoulombop_" +
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

  /// emits the element-wise irrep combination of the children, e.g.
  /// target = 2.0*src0 - src1 - src2 for SymTLDiagB
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    rkb_fold_generate_code_from_terms(
        *this, fold_terms(target_->oper()->descr().component_index()),
        "CR_11_σpCoulombσp_11", context, dims, funcname, decl, def);
  }
};

template <typename F>
CR_11_σpCoulombσp_11<F>::CR_11_σpCoulombσp_11(
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

}  // CR_11_σpCoulombσp_11<F>::CR_11_σpCoulombσp_11

}  // namespace libint2

#endif  // LIBINT_COMP_11_ΣPCOULOMBΣP_11_H
