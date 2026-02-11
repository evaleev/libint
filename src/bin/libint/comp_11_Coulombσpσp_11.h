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

#include <gaussoper.h>
#include <generic_rr.h>
#include <twoprep_11_11.h>

namespace libint2 {

/**
 * this computes integral of
 * \f$ \frac{1}{r_{ij}} \sigma \cdot \hat{p}_1 \sigma \cdot \hat{p}_2 \f$ over
 * CGShell/CGF by rewriting it as a linear combination of integrals over
 * derivatives of \frac{1}{r_{ij}}
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
  CR_11_Coulombσpσp_11(const std::shared_ptr<TargetType> &, unsigned int = 0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_11_Coulombσpσp_11<F>::CR_11_Coulombσpσp_11(
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

  // can express integrals of Coulombσpσp in terms of
  // derivative integrals of 1/r12 for primitive Gaussians
  // only
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

  F c_x{c};
  c_x.deriv().inc(x);  // d(c)/dx = c_x
  F c_y{c};
  c_y.deriv().inc(y);  // d(c)/dy = c_y
  F c_z{c};
  c_z.deriv().inc(z);  // d(c)/dz = c_z

  F d_x{d};
  d_x.deriv().inc(x);  // d(d)/dx = d_x
  F d_y{d};
  d_y.deriv().inc(y);  // d(d)/dy = d_y
  F d_z{d};
  d_z.deriv().inc(z);  // d(d)/dz = d_z

  // Component wise generation for quaternion ( a b | 1/r12 | (σ.p) c (σ.p) d )
  switch (oper->descr().quaternion_index()) {
    case 0: {
      // zeroth component = (a b | c_x d_x) + (a b | c_y d_y) + (a b | c_z d_z)
      auto a_b_cx_dx = factory.make_child(a, b, c_x, d_x, zero_m);
      auto a_b_cy_dy = factory.make_child(a, b, c_y, d_y, zero_m);
      auto a_b_cz_dz = factory.make_child(a, b, c_z, d_z, zero_m);
      if (is_simple()) {
        expr_ = a_b_cx_dx + a_b_cy_dy + a_b_cz_dz;
        nflops_ += 2;
      }
    } break;
    case 1: {
      // x component = (a b | c_y d_z) - (a b | c_z d_y)
      auto a_b_cy_dz = factory.make_child(a, b, c_y, d_z, zero_m);
      auto a_b_cz_dy = factory.make_child(a, b, c_z, d_y, zero_m);
      if (is_simple()) {
        expr_ = a_b_cy_dz - a_b_cz_dy;
        nflops_ += 1;
      }
    } break;
    case 2: {
      // y component = (a b | c_z d_x) - (a b | c_x d_z)
      auto a_b_cz_dx = factory.make_child(a, b, c_z, d_x, zero_m);
      auto a_b_cx_dz = factory.make_child(a, b, c_x, d_z, zero_m);
      if (is_simple()) {
        expr_ = a_b_cz_dx - a_b_cx_dz;
        nflops_ += 1;
      }
    } break;
    case 3: {
      // z component = (a b | c_x d_y) - (a b | c_y d_x)
      auto a_b_cx_dy = factory.make_child(a, b, c_x, d_y, zero_m);
      auto a_b_cy_dx = factory.make_child(a, b, c_y, d_x, zero_m);
      if (is_simple()) {
        expr_ = a_b_cx_dy - a_b_cy_dx;
        nflops_ += 1;
      }
    } break;
    default:
      throw std::runtime_error(
          "CR_11_Coulombσpσp_11: invalid quaternionic index");
  }

}  // CR_11_Coulombσpσp_11<F>::CR_11_Coulombσpσp_11
};  // namespace libint2

#endif  // LIBINT_COMP_11_COULOMBΣPΣP_11_H
