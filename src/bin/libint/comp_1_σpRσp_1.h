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

#ifndef LIBINT_COMP_1_ΣPRΣP_1_H
#define LIBINT_COMP_1_ΣPRΣP_1_H

#include <generic_rr.h>

namespace libint2 {

/**
 * Computes the integral of \f$ \sigma \cdot \hat{p}\, r_k\, \sigma \cdot
 * \hat{p} \f$ over CGShell/CGF by folding the 9 raw \f$ \sigma_a \partial_a r_k
 * \sigma_b \partial_b \f$ dyadics per dipole direction \f$ k \f$ to 4
 * Pauli-quaternion components via \f$ \sigma_a \sigma_b = \delta_{ab} +
 * i\epsilon_{abc}\sigma_c \f$. 12 outputs total = 3 dipole directions × 4
 * Pauli components, mirroring the σpVσp fold but with the central operator
 * being a Cartesian dipole instead of the electrostatic potential V.
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_1_σpRσp_1
    : public GenericRecurrenceRelation<
          CR_1_σpRσp_1<F>, F, GenIntegralSet_1_1<F, σpRσpOper, EmptySet>> {
 public:
  typedef CR_1_σpRσp_1<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpRσpOper OperType;
  typedef GenIntegralSet_1_1<F, σpRσpOper, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 100;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  CR_1_σpRσp_1(const std::shared_ptr<TargetType> &, unsigned int = 0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_1_σpRσp_1<F>::CR_1_σpRσp_1(const std::shared_ptr<TargetType> &Tint,
                              unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  const auto &a = Tint->bra(0, 0);
  const auto &b = Tint->ket(0, 0);
  const auto &oper = Tint->oper();

  // express σ·p r_k σ·p in terms of derivative integrals of the dipole
  // operator r_k for primitive Gaussians only
  if (a.contracted() || b.contracted()) return;

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  ChildFactory<ThisType,
               GenIntegralSet_1_1<BasisFunctionType, CartesianMultipoleOper<3u>,
                                  EmptySet>>
      factory(this);

  constexpr auto x = 0;
  constexpr auto y = 1;
  constexpr auto z = 2;

  F Dx_a{a};
  Dx_a.deriv().inc(x);
  F Dx_b{b};
  Dx_b.deriv().inc(x);
  F Dy_a{a};
  Dy_a.deriv().inc(y);
  F Dy_b{b};
  Dy_b.deriv().inc(y);
  F Dz_a{a};
  Dz_a.deriv().inc(z);
  F Dz_b{b};
  Dz_b.deriv().inc(z);

  // Build the dipole multipole descriptor for direction k.
  const auto k = oper->descr().dipole_dir();
  CartesianMultipole_Descr<3u> mu_k;
  mu_k.inc(k, 1);  // r_k = (kx,ky,kz) with k_k = 1, others 0

  // Pauli quaternion fold per (k, q):
  //   q=0: trace δ_ab → Σ_a (∂_a μ | r_k | ∂_a ν)
  //   q=1: σ_x antisym → (∂_y μ | r_k | ∂_z ν) − (∂_z μ | r_k | ∂_y ν)
  //   q=2: σ_y antisym → (∂_z μ | r_k | ∂_x ν) − (∂_x μ | r_k | ∂_z ν)
  //   q=3: σ_z antisym → (∂_x μ | r_k | ∂_y ν) − (∂_y μ | r_k | ∂_x ν)
  switch (oper->descr().quaternion_index()) {
    case 0: {
      auto Dx_a_R_Dx_b = factory.make_child(Dx_a, Dx_b, EmptySet(), mu_k);
      auto Dy_a_R_Dy_b = factory.make_child(Dy_a, Dy_b, EmptySet(), mu_k);
      auto Dz_a_R_Dz_b = factory.make_child(Dz_a, Dz_b, EmptySet(), mu_k);
      if (is_simple()) {
        expr_ = Dx_a_R_Dx_b + Dy_a_R_Dy_b + Dz_a_R_Dz_b;
        nflops_ += 2;
      }
    } break;
    case 1: {
      auto Dy_a_R_Dz_b = factory.make_child(Dy_a, Dz_b, EmptySet(), mu_k);
      auto Dz_a_R_Dy_b = factory.make_child(Dz_a, Dy_b, EmptySet(), mu_k);
      if (is_simple()) {
        expr_ = Dy_a_R_Dz_b - Dz_a_R_Dy_b;
        nflops_ += 1;
      }
    } break;
    case 2: {
      auto Dz_a_R_Dx_b = factory.make_child(Dz_a, Dx_b, EmptySet(), mu_k);
      auto Dx_a_R_Dz_b = factory.make_child(Dx_a, Dz_b, EmptySet(), mu_k);
      if (is_simple()) {
        expr_ = Dz_a_R_Dx_b - Dx_a_R_Dz_b;
        nflops_ += 1;
      }
    } break;
    case 3: {
      auto Dx_a_R_Dy_b = factory.make_child(Dx_a, Dy_b, EmptySet(), mu_k);
      auto Dy_a_R_Dx_b = factory.make_child(Dy_a, Dx_b, EmptySet(), mu_k);
      if (is_simple()) {
        expr_ = Dx_a_R_Dy_b - Dy_a_R_Dx_b;
        nflops_ += 1;
      }
    } break;
    default:
      throw std::runtime_error("CR_1_σpRσp_1: invalid quaternionic index");
  }

}  // CR_1_σpRσp_1<F>::CR_1_σpRσp_1

};  // namespace libint2

#endif  // LIBINT_COMP_1_ΣPRΣP_1_H
