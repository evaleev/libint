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

#ifndef LIBINT_COMP_1_ΣPVΣP_1_H
#define LIBINT_COMP_1_ΣPVΣP_1_H

#include <generic_rr.h>
#include <rkb_fold_1body.h>

namespace libint2 {

/**
 * this computes integral of
 * \f$ \sigma \cdot \hat{p} V \sigma \cdot \hat{p} \f$ over CGShell/CGF
 * by rewriting it as a linear combination of integrals over electrostatic
 * potential \f$ V \f$
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_1_σpVσp_1
    : public GenericRecurrenceRelation<
          CR_1_σpVσp_1<F>, F, GenIntegralSet_1_1<F, σpVσpOper, EmptySet>> {
 public:
  typedef CR_1_σpVσp_1<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpVσpOper OperType;
  typedef GenIntegralSet_1_1<F, σpVσpOper, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 100;  // TODO figure out

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_1_σpVσp_1(const std::shared_ptr<TargetType> &, unsigned int = 0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_1_σpVσp_1<F>::CR_1_σpVσp_1(const std::shared_ptr<TargetType> &Tint,
                              unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  const auto &a = Tint->bra(0, 0);
  const auto &b = Tint->ket(0, 0);
  const auto &oper = Tint->oper();

  // can express integrals of σpVσp in terms of derivative integrals of V for
  // primitive Gaussians only
  if (a.contracted() || b.contracted()) return;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_1_1<BasisFunctionType, ElecPotOper, mType>>
      factory(this);

  // σ·p V σ·p fold: inner child is the electrostatic-potential integral
  // (Da | V | Db). Shared structure lives in sigma_p_1body_fold().
  sigma_p_1body_fold(
      a, b, oper->descr().quaternion_index(), is_simple(), expr_, nflops_,
      [&](const F &Da, const F &Db) {
        return factory.make_child(Da, Db, zero_m);
      });

}  // CR_1_σpVσp_1<F>::CR_1_σpVσp_1

};  // namespace libint2

#endif  // LIBINT_COMP_1_ΣPVΣP_1_H
