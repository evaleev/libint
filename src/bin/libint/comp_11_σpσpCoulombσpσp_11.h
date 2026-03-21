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
 * Computes integral of
 * \f$ (\sigma_1 \cdot \hat{p}_a)(\sigma_1 \cdot \hat{p}_b)
 *     \frac{1}{r_{12}}
 *     (\sigma_2 \cdot \hat{p}_c)(\sigma_2 \cdot \hat{p}_d) \f$
 * over CGShell/CGF by rewriting it as a linear combination of integrals
 * over derivatives of \f$ \frac{1}{r_{12}} \f$.
 *
 * The two sigma operators act on independent spin spaces (electron 1 and
 * electron 2). Using the Dirac identity (see e.g. Eq. 1.27 of I. P. Grant,
 * "Relativistic Quantum Theory of Atoms and Molecules", Springer, 2007):
 *   \f$ (\sigma \cdot a)(\sigma \cdot b) = (a \cdot b)I
 *       + i\sigma \cdot (a \times b) \f$
 * applied independently to each particle's spin space gives a tensor product
 * of two quaternions with \f$ 4 \times 4 = 16 \f$ components:
 *
 *   index = 4 * bra_spin_index + ket_spin_index
 *
 * where spin indices are: 0=S (scalar/dot product), 1=X, 2=Y, 3=Z
 * (cross product components).
 *
 * The 16 components map to:
 *   T1 (index 0):      SS = (a.b)(c.d)                  [scalar x scalar]
 *   T2 (indices 1-3):  SX,SY,SZ = (a.b)(cxd)_{x,y,z}   [scalar x spin]
 *   T3 (indices 4-6):  XS,YS,ZS = (axb)_{x,y,z}(c.d)   [spin x scalar]
 *   T4 (indices 7-15): XX..ZZ = -(axb)_i(cxd)_j         [spin x spin]
 *
 * Sign convention: T4 components include the minus sign from \f$ i^2 = -1 \f$,
 * arising from the product of two \f$ i \f$ factors in the Dirac identity:
 *   \f$ [i\sigma_1 \cdot (a \times b)] \otimes [i\sigma_2 \cdot (c \times d)]
 *     = -\sigma_{1,i} \otimes \sigma_{2,j}\; (a \times b)_i\, (c \times d)_j
 *   \f$
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

  // 16-component generation for two independent spin spaces:
  // ( (σ₁.p) a (σ₁.p) b | 1/r12 | (σ₂.p) c (σ₂.p) d )
  //
  // Option A (tensor product) ordering: index = 4 * bra_spin + ket_spin
  //   bra_spin = index / 4,  ket_spin = index % 4
  //   spin indices: S=0, X=1, Y=2, Z=3
  //
  // Row bra=S: 0=SS, 1=SX, 2=SY, 3=SZ          (T1 + T2)
  // Row bra=X: 4=XS, 5=XX, 6=XY, 7=XZ          (T3 + T4)
  // Row bra=Y: 8=YS, 9=YX, 10=YY, 11=YZ        (T3 + T4)
  // Row bra=Z: 12=ZS, 13=ZX, 14=ZY, 15=ZZ      (T3 + T4)
  //
  // T4 components include minus sign from i^2 = -1.
  switch (oper->descr().quaternion_index()) {
    // ===== 0: SS = (a.b)(c.d) =====
    case 0: {
      auto xxxx = mc(x, x, x, x);
      auto xxyy = mc(x, x, y, y);
      auto xxzz = mc(x, x, z, z);
      auto yyxx = mc(y, y, x, x);
      auto yyyy = mc(y, y, y, y);
      auto yyzz = mc(y, y, z, z);
      auto zzxx = mc(z, z, x, x);
      auto zzyy = mc(z, z, y, y);
      auto zzzz = mc(z, z, z, z);
      if (is_simple()) {
        expr_ = xxxx + xxyy + xxzz + yyxx + yyyy + yyzz + zzxx + zzyy + zzzz;
        nflops_ += 8;
      }
    } break;
    // ===== 1: SX = (a.b)(c×d)_x =====
    case 1: {
      auto xxyz = mc(x, x, y, z);
      auto xxzy = mc(x, x, z, y);
      auto yyyz = mc(y, y, y, z);
      auto yyzy = mc(y, y, z, y);
      auto zzyz = mc(z, z, y, z);
      auto zzzy = mc(z, z, z, y);
      if (is_simple()) {
        expr_ = xxyz - xxzy + yyyz - yyzy + zzyz - zzzy;
        nflops_ += 5;
      }
    } break;
    // ===== 2: SY = (a.b)(c×d)_y =====
    case 2: {
      auto xxzx = mc(x, x, z, x);
      auto xxxz = mc(x, x, x, z);
      auto yyzx = mc(y, y, z, x);
      auto yyxz = mc(y, y, x, z);
      auto zzzx = mc(z, z, z, x);
      auto zzxz = mc(z, z, x, z);
      if (is_simple()) {
        expr_ = xxzx - xxxz + yyzx - yyxz + zzzx - zzxz;
        nflops_ += 5;
      }
    } break;
    // ===== 3: SZ = (a.b)(c×d)_z =====
    case 3: {
      auto xxxy = mc(x, x, x, y);
      auto xxyx = mc(x, x, y, x);
      auto yyxy = mc(y, y, x, y);
      auto yyyx = mc(y, y, y, x);
      auto zzxy = mc(z, z, x, y);
      auto zzyx = mc(z, z, y, x);
      if (is_simple()) {
        expr_ = xxxy - xxyx + yyxy - yyyx + zzxy - zzyx;
        nflops_ += 5;
      }
    } break;
    // ===== 4: XS = (a×b)_x(c.d) =====
    case 4: {
      auto yzxx = mc(y, z, x, x);
      auto zyxx = mc(z, y, x, x);
      auto yzyy = mc(y, z, y, y);
      auto zyyy = mc(z, y, y, y);
      auto yzzz = mc(y, z, z, z);
      auto zyzz = mc(z, y, z, z);
      if (is_simple()) {
        expr_ = yzxx - zyxx + yzyy - zyyy + yzzz - zyzz;
        nflops_ += 5;
      }
    } break;
    // ===== 5: XX = -(a×b)_x(c×d)_x (minus from i²=-1) =====
    case 5: {
      auto yzyz = mc(y, z, y, z);
      auto yzzy = mc(y, z, z, y);
      auto zyyz = mc(z, y, y, z);
      auto zyzy = mc(z, y, z, y);
      if (is_simple()) {
        expr_ = yzzy - yzyz + zyyz - zyzy;
        nflops_ += 3;
      }
    } break;
    // ===== 6: XY = -(a×b)_x(c×d)_y =====
    case 6: {
      auto yzzx = mc(y, z, z, x);
      auto yzxz = mc(y, z, x, z);
      auto zyzx = mc(z, y, z, x);
      auto zyxz = mc(z, y, x, z);
      if (is_simple()) {
        expr_ = yzxz - yzzx + zyzx - zyxz;
        nflops_ += 3;
      }
    } break;
    // ===== 7: XZ = -(a×b)_x(c×d)_z =====
    case 7: {
      auto yzxy = mc(y, z, x, y);
      auto yzyx = mc(y, z, y, x);
      auto zyxy = mc(z, y, x, y);
      auto zyyx = mc(z, y, y, x);
      if (is_simple()) {
        expr_ = yzyx - yzxy + zyxy - zyyx;
        nflops_ += 3;
      }
    } break;
    // ===== 8: YS = (a×b)_y(c.d) =====
    case 8: {
      auto zxxx = mc(z, x, x, x);
      auto xzxx = mc(x, z, x, x);
      auto zxyy = mc(z, x, y, y);
      auto xzyy = mc(x, z, y, y);
      auto zxzz = mc(z, x, z, z);
      auto xzzz = mc(x, z, z, z);
      if (is_simple()) {
        expr_ = zxxx - xzxx + zxyy - xzyy + zxzz - xzzz;
        nflops_ += 5;
      }
    } break;
    // ===== 9: YX = -(a×b)_y(c×d)_x =====
    case 9: {
      auto zxyz = mc(z, x, y, z);
      auto zxzy = mc(z, x, z, y);
      auto xzyz = mc(x, z, y, z);
      auto xzzy = mc(x, z, z, y);
      if (is_simple()) {
        expr_ = zxzy - zxyz + xzyz - xzzy;
        nflops_ += 3;
      }
    } break;
    // ===== 10: YY = -(a×b)_y(c×d)_y =====
    case 10: {
      auto zxzx = mc(z, x, z, x);
      auto zxxz = mc(z, x, x, z);
      auto xzzx = mc(x, z, z, x);
      auto xzxz = mc(x, z, x, z);
      if (is_simple()) {
        expr_ = zxxz - zxzx + xzzx - xzxz;
        nflops_ += 3;
      }
    } break;
    // ===== 11: YZ = -(a×b)_y(c×d)_z =====
    case 11: {
      auto zxxy = mc(z, x, x, y);
      auto zxyx = mc(z, x, y, x);
      auto xzxy = mc(x, z, x, y);
      auto xzyx = mc(x, z, y, x);
      if (is_simple()) {
        expr_ = zxyx - zxxy + xzxy - xzyx;
        nflops_ += 3;
      }
    } break;
    // ===== 12: ZS = (a×b)_z(c.d) =====
    case 12: {
      auto xyxx = mc(x, y, x, x);
      auto yxxx = mc(y, x, x, x);
      auto xyyy = mc(x, y, y, y);
      auto yxyy = mc(y, x, y, y);
      auto xyzz = mc(x, y, z, z);
      auto yxzz = mc(y, x, z, z);
      if (is_simple()) {
        expr_ = xyxx - yxxx + xyyy - yxyy + xyzz - yxzz;
        nflops_ += 5;
      }
    } break;
    // ===== 13: ZX = -(a×b)_z(c×d)_x =====
    case 13: {
      auto xyyz = mc(x, y, y, z);
      auto xyzy = mc(x, y, z, y);
      auto yxyz = mc(y, x, y, z);
      auto yxzy = mc(y, x, z, y);
      if (is_simple()) {
        expr_ = xyzy - xyyz + yxyz - yxzy;
        nflops_ += 3;
      }
    } break;
    // ===== 14: ZY = -(a×b)_z(c×d)_y =====
    case 14: {
      auto xyzx = mc(x, y, z, x);
      auto xyxz = mc(x, y, x, z);
      auto yxzx = mc(y, x, z, x);
      auto yxxz = mc(y, x, x, z);
      if (is_simple()) {
        expr_ = xyxz - xyzx + yxzx - yxxz;
        nflops_ += 3;
      }
    } break;
    // ===== 15: ZZ = -(a×b)_z(c×d)_z =====
    case 15: {
      auto xyxy = mc(x, y, x, y);
      auto xyyx = mc(x, y, y, x);
      auto yxxy = mc(y, x, x, y);
      auto yxyx = mc(y, x, y, x);
      if (is_simple()) {
        expr_ = xyyx - xyxy + yxxy - yxyx;
        nflops_ += 3;
      }
    } break;
    default:
      throw std::runtime_error(
          "CR_11_σpσpCoulombσpσp_11: invalid component index (expected 0-15)");
  }

}  // CR_11_σpσpCoulombσpσp_11<F>::CR_11_σpσpCoulombσpσp_11
};  // namespace libint2

#endif  // LIBINT_COMP_11_ΣPΣPCOULOMBΣPΣP_11_H
