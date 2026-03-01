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

#include <gaussoper.h>
#include <generic_rr.h>
#include <twoprep_11_11.h>

namespace libint2 {

/**
 * this computes integral of
 * \sigma \cdot \hat{p}_1 \sigma \cdot \hat{p}_2 \f$ \frac{1}{r_{ij}} \sigma
 * \cdot \hat{p}_3 \sigma \cdot \hat{p}_4 \f$ over CGShell/CGF by rewriting it
 * as a linear combination of integrals over derivatives of \frac{1}{r_{ij}}
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
  CR_11_σpσpCoulombσpσp_11(const std::shared_ptr<TargetType> &,
                           unsigned int = 0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_11_σpσpCoulombσpσp_11<F>::CR_11_σpσpCoulombσpσp_11(
    const std::shared_ptr<TargetType> &Tint, unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_bra(/* particle */ 1) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 1) == 1);

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  const auto &oper = Tint->oper();

  if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
    return;

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_11_11<BasisFunctionType, TwoPRep, mType>>
      factory(this);

  constexpr auto x = 0;
  constexpr auto y = 1;
  constexpr auto z = 2;

  auto mc = [&](const int r1, const int r2, const int r3, const int r4) {
    F a_r1{a};
    a_r1.deriv().inc(r1);
    F b_r2{b};
    b_r2.deriv().inc(r2);
    F c_r3{c};
    c_r3.deriv().inc(r3);
    F d_r4{d};
    d_r4.deriv().inc(r4);
    return factory.make_child(a_r1, b_r2, c_r3, d_r4, zero_m);
  };

  // Component wise generation for quaternion :
  // ( (σ.p) a (σ.p) b | 1/r12 | (σ.p) c (σ.p) d )
  switch (oper->descr().quaternion_index()) {
    case 0: {
      auto xxxx = mc(x, x, x, x);
      auto yyxx = mc(y, y, x, x);
      auto zzxx = mc(z, z, x, x);
      auto yxyx = mc(y, x, y, x);
      auto xyyx = mc(x, y, y, x);
      auto yxxy = mc(y, x, x, y);
      auto xyxy = mc(x, y, x, y);
      auto xxyy = mc(x, x, y, y);
      auto yyyy = mc(y, y, y, y);
      auto zzyy = mc(z, z, y, y);
      auto xxzz = mc(x, x, z, z);
      auto yyzz = mc(y, y, z, z);
      auto zzzz = mc(z, z, z, z);
      if (is_simple()) {
        expr_ = xxxx + yyxx + zzxx - yxyx + xyyx + yxxy - xyxy + xxyy + yyyy +
                zzyy + xxzz + yyzz + zzzz;
        nflops_ += 12;
      }
    } break;
    case 1: {
      auto zxzx = mc(z, x, z, x);
      auto xzzx = mc(x, z, z, x);
      auto zyzy = mc(z, y, z, y);
      auto yzzy = mc(y, z, z, y);
      auto zxxz = mc(z, x, x, z);
      auto xzxz = mc(x, z, x, z);
      auto zyyz = mc(z, y, y, z);
      auto yzyz = mc(y, z, y, z);
      if (is_simple()) {
        expr_ = zxzx - xzzx - zyzy + yzzy - zxxz + xzxz + zyyz - yzyz;
        nflops_ += 7;
      }
    } break;
    case 2: {
      auto zyzx = mc(z, y, z, x);
      auto yzzx = mc(y, z, z, x);
      auto zxzy = mc(z, x, z, y);
      auto xzzy = mc(x, z, z, y);
      auto zyxz = mc(z, y, x, z);
      auto yzxz = mc(y, z, x, z);
      auto zxyz = mc(z, x, y, z);
      auto xzyz = mc(x, z, y, z);
      if (is_simple()) {
        // swapped order of first two terms compiler does not like negative sign
        // in front of first term
        expr_ = yzzx - zyzx - zxzy + xzzy + zyxz - yzxz + zxyz - xzyz;
        nflops_ += 7;
      }
    } break;
    case 3: {
      auto yxxx = mc(y, x, x, x);
      auto xyxx = mc(x, y, x, x);
      auto xxyx = mc(x, x, y, x);
      auto yyyx = mc(y, y, y, x);
      auto zzyx = mc(z, z, y, x);
      auto xxxy = mc(x, x, x, y);
      auto yyxy = mc(y, y, x, y);
      auto zzxy = mc(z, z, x, y);
      auto yxyy = mc(y, x, y, y);
      auto xyyy = mc(x, y, y, y);
      auto yxzz = mc(y, x, z, z);
      auto xyzz = mc(x, y, z, z);
      if (is_simple()) {
        expr_ = xyxx - yxxx - xxyx - yyyx - zzyx + xxxy + yyxy + zzxy - yxyy +
                xyyy - yxzz + xyzz;
        nflops_ += 11;
      }
    } break;
    default:
      throw std::runtime_error(
          "CR_11_σpσpCoulombσpσp_11: invalid quaternionic index");
  }

}  // CR_11_σpσpCoulombσpσp_11<F>::CR_11_σpσpCoulombσpσp_11
};  // namespace libint2

#endif  // LIBINT_COMP_11_ΣPΣPCOULOMBΣPΣP_11_H
