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
#include <rkb_fold_1body.h>

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
class CR_1_σpeμσp_1
    : public GenericRecurrenceRelation<
          CR_1_σpeμσp_1<F>, F, GenIntegralSet_1_1<F, σpeμσpOper, EmptySet>> {
 public:
  typedef CR_1_σpeμσp_1<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpeμσpOper OperType;
  typedef GenIntegralSet_1_1<F, σpeμσpOper, EmptySet> TargetType;
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

  CR_1_σpeμσp_1(const std::shared_ptr<TargetType> &, unsigned int = 0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_1_σpeμσp_1<F>::CR_1_σpeμσp_1(const std::shared_ptr<TargetType> &Tint,
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

  ChildFactory<ThisType,
               GenIntegralSet_1_1<BasisFunctionType, CartesianMultipoleOper<3u>,
                                  EmptySet>>
      factory(this);

  // Build the dipole multipole descriptor for direction k.
  const auto k = oper->descr().dipole_dir();
  CartesianMultipole_Descr<3u> mu_k;
  mu_k.inc(k, 1);  // r_k = (kx,ky,kz) with k_k = 1, others 0

  // σ·p r_k σ·p fold: inner child is the Cartesian-dipole integral
  // (Da | r_k | Db). Shared structure lives in sigma_p_1body_fold().
  sigma_p_1body_fold(
      a, b, oper->descr().quaternion_index(), is_simple(), expr_, nflops_,
      [&](const F &Da, const F &Db) {
        return factory.make_child(Da, Db, EmptySet(), mu_k);
      });

}  // CR_1_σpeμσp_1<F>::CR_1_σpeμσp_1

};  // namespace libint2

#endif  // LIBINT_COMP_1_ΣPRΣP_1_H
