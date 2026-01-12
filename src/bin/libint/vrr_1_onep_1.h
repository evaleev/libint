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

#ifndef _libint2_src_bin_libint_vrr1onep1_h_
#define _libint2_src_bin_libint_vrr1onep1_h_

#include <generic_rr.h>
#include <onep_1_1.h>

#include <cassert>

namespace libint2 {

/** VRR Recurrence Relation for 1-e overlap integrals.
 * \tparam where specifies whether angular momentum is decreased, in bra or ket.
 */
template <class BFSet, FunctionPosition where>
class VRR_1_Overlap_1 : public GenericRecurrenceRelation<
                            VRR_1_Overlap_1<BFSet, where>, BFSet,
                            GenIntegralSet_1_1<BFSet, OverlapOper, EmptySet> > {
 public:
  typedef VRR_1_Overlap_1 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_1_1<BFSet, OverlapOper, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 9;

  using ParentType::Instance;

  /// Default directionality
  static bool directional() { return ParentType::default_directional(); }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  VRR_1_Overlap_1(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "OSVRROverlap"; }
};

template <class F, FunctionPosition where>
VRR_1_Overlap_1<F, where>::VRR_1_Overlap_1(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  const F& _1 = unit<F>(dir);

  {  // can't apply to contracted basis functions
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    if (a.contracted() || b.contracted()) return;
  }

  // if derivative integrals, there will be extra terms (Eq. (143) in Obara &
  // Saika JCP 89)
  const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero();

  typedef TargetType ChildType;
  ChildFactory<ThisType, ChildType> factory(this);

  // Build on A or B
  {
    // bf quantum on the build center subtracted by 1
    auto a = (where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    if (!exists(a)) return;
    auto b = (where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
    if (!exists(b)) return;

    auto AB = factory.make_child(a, b);
    if (is_simple()) {
      expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB;
      nflops_ += 1;
    }

    auto am1 = a - _1;
    auto am1_exists = exists(am1);
    auto bm1 = b - _1;
    auto bm1_exists = exists(bm1);

    if (am1_exists) {
      auto Am1B = factory.make_child(am1, b);
      if (is_simple()) {
        expr_ += (Scalar(a[dir]) * Scalar("oo2z")) * Am1B;
        nflops_ += 3;
      }
    }
    if (bm1_exists) {
      auto ABm1 = factory.make_child(a, bm1);
      if (is_simple()) {
        expr_ += (Scalar(b[dir]) * Scalar("oo2z")) * ABm1;
        nflops_ += 3;
      }
    }
  }

  // if got here, can decrement by at least 1 quantum
  // add additional derivative terms
  if (deriv) {
    // bf quantum on the build center subtracted by 1
    F a(where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    F b(where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));

    // treatment of derivative terms differs for shell sets and integrals
    // since in computing shell sets transfer/build will occur in all 3
    // directions change in up to all three derivative indices will occur
    for (unsigned int dxyz = 0; dxyz < 3; ++dxyz) {
      if (is_simple() && dxyz != dir)  // for integrals only consider
                                       // derivatives in THE build direction
        continue;

      OriginDerivative<3u> _d1;
      _d1.inc(dxyz);

      std::shared_ptr<DGVertex> _nullptr;

      // dA - _1?
      {
        const OriginDerivative<3u> dAm1(dA - _d1);
        if (exists(dAm1)) {  // yes
          a.deriv() = dAm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ -= Vector(dA)[dxyz] * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ += Vector(dA)[dxyz] * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          a.deriv() = dA;
        }
      }

      // dB - _1?
      {
        const OriginDerivative<3u> dBm1(dB - _d1);
        if (exists(dBm1)) {  // yes
          b.deriv() = dBm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ += Vector(dB)[dxyz] * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ -= Vector(dB)[dxyz] * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          b.deriv() = dB;
        }
      }
    }
  }  // end of deriv

  return;
}

/** VRR Recurrence Relation for 1-d overlap integrals.
 * \tparam where specifies whether quantum number is decreased, in bra or ket.
 */
template <CartesianAxis Axis, FunctionPosition where>
class VRR_1_Overlap_1_1d
    : public GenericRecurrenceRelation<
          VRR_1_Overlap_1_1d<Axis, where>, CGF1d<Axis>,
          GenIntegralSet_1_1<CGF1d<Axis>, OverlapOper, EmptySet> > {
 public:
  typedef VRR_1_Overlap_1_1d ThisType;
  typedef CGF1d<Axis> BasisFunctionType;
  typedef GenIntegralSet_1_1<BasisFunctionType, OverlapOper, EmptySet>
      TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 9;
  static constexpr auto axis = Axis;

  using ParentType::Instance;

  /// Default directionality
  static bool directional() { return ParentType::default_directional(); }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  VRR_1_Overlap_1_1d(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() {
    return std::string("OSVRROverlap") + to_string(axis);
  }
};

template <CartesianAxis Axis, FunctionPosition where>
VRR_1_Overlap_1_1d<Axis, where>::VRR_1_Overlap_1_1d(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  assert(dir == 0);  // this integral is along 1 axis only

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  typedef CGF1d<Axis> F;
  const F& _1 = unit<F>(dir);

  {  // can't apply to contracted basis functions
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    if (a.contracted() || b.contracted()) return;
  }

  // if derivative integrals, there will be extra terms (Eq. (143) in Obara &
  // Saika JCP 89)
  const OriginDerivative<1u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<1u> dB = Tint->ket(0, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero();

  typedef TargetType ChildType;
  ChildFactory<ThisType, ChildType> factory(this);

  // handle the special case of (0|0) integral
  // to avoid complications with non-precomputed shell blocks it's "computed
  // by copying from inteval
  auto zero = Tint->bra(0, 0).zero() and Tint->ket(0, 0).zero() and not deriv;
  if (zero) {
    std::shared_ptr<DGVertex> int00 = Vector("_0_Overlap_0")[Axis];
    expr_ = Scalar(0u) + int00;
    this->add_child(int00);
    return;
  }

  // Build on A or B
  {
    // bf quantum on the build center subtracted by 1
    auto a = (where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    if (!exists(a)) return;
    auto b = (where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
    if (!exists(b)) return;

    auto AB = factory.make_child(a, b);
    if (is_simple()) {
      expr_ = Vector(where == InBra ? "PA" : "PB")[Axis] * AB;
      nflops_ += 1;
    }

    auto am1 = a - _1;
    auto am1_exists = exists(am1);
    auto bm1 = b - _1;
    auto bm1_exists = exists(bm1);

    if (am1_exists) {
      auto Am1B = factory.make_child(am1, b);
      if (is_simple()) {
        expr_ += (Scalar(a.qn()) * Scalar("oo2z")) * Am1B;
        nflops_ += 3;
      }
    }
    if (bm1_exists) {
      auto ABm1 = factory.make_child(a, bm1);
      if (is_simple()) {
        expr_ += (Scalar(b.qn()) * Scalar("oo2z")) * ABm1;
        nflops_ += 3;
      }
    }
  }

  // if got here, can decrement by at least 1 quantum
  // add additional derivative terms
  if (deriv) {
    // bf quantum on the build center subtracted by 1
    F a(where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    F b(where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));

    {
      OriginDerivative<1u> _d1;
      _d1.inc(0);

      std::shared_ptr<DGVertex> _nullptr;

      // dA - _1?
      {
        const OriginDerivative<1u> dAm1(dA - _d1);
        if (exists(dAm1)) {  // yes
          a.deriv() = dAm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ -= Scalar(dA[0]) * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ += Scalar(dA[0]) * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          a.deriv() = dA;
        }
      }

      // dB - _1?
      {
        const OriginDerivative<1u> dBm1(dB - _d1);
        if (exists(dBm1)) {  // yes
          b.deriv() = dBm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ += Scalar(dB[0]) * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ -= Scalar(dB[0]) * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          b.deriv() = dB;
        }
      }
    }
  }  // end of deriv

  return;
}

/** VRR Recurrence Relation for 1-e kinetic energy integrals.
 * \tparam where specifies whether angular momentum is decreased, in bra or ket.
 */
template <class BFSet, FunctionPosition where>
class VRR_1_Kinetic_1 : public GenericRecurrenceRelation<
                            VRR_1_Kinetic_1<BFSet, where>, BFSet,
                            GenIntegralSet_1_1<BFSet, KineticOper, EmptySet> > {
 public:
  typedef VRR_1_Kinetic_1 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_1_1<BFSet, KineticOper, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 9;

  using ParentType::Instance;

  /// Default directionality
  static bool directional() { return ParentType::default_directional(); }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  VRR_1_Kinetic_1(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "OSVRRKinetic"; }
};

template <class F, FunctionPosition where>
VRR_1_Kinetic_1<F, where>::VRR_1_Kinetic_1(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  const F& _1 = unit<F>(dir);

  {  // can't apply to contracted basis functions
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    if (a.contracted() || b.contracted()) return;
  }

  // if derivative integrals, there will be extra terms (Eq. (143) in Obara &
  // Saika JCP 89)
  const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero();

  typedef TargetType Child1Type;
  ChildFactory<ThisType, Child1Type> factory(this);
  typedef GenIntegralSet_1_1<F, OverlapOper, EmptySet> Child2Type;
  ChildFactory<ThisType, Child2Type> overlap_factory(this);

  // Build on A or B
  {
    // bf quantum on the build center subtracted by 1
    auto a = (where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    if (!exists(a)) return;
    auto b = (where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
    if (!exists(b)) return;

    auto AB = factory.make_child(a, b);
    if (is_simple()) {
      expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB;
      nflops_ += 1;
    }

    auto am1 = a - _1;
    if (exists(am1)) {
      auto Am1B = factory.make_child(am1, b);
      auto S_Am1B = (where == InBra) ? overlap_factory.make_child(am1, b)
                                     : std::shared_ptr<DGVertex>();
      if (is_simple()) {
        if (where == InBra) {
          expr_ += Scalar(a[dir]) * (Scalar("oo2z") * Am1B -
                                     Scalar("rho12_over_alpha1") * S_Am1B);
          nflops_ += 5;
        } else {
          expr_ += Scalar(a[dir]) * Scalar("oo2z") * Am1B;
          nflops_ += 3;
        }
      }
    }
    auto bm1 = b - _1;
    if (exists(bm1)) {
      auto ABm1 = factory.make_child(a, bm1);
      auto S_ABm1 = (where == InKet) ? overlap_factory.make_child(a, bm1)
                                     : std::shared_ptr<DGVertex>();
      if (is_simple()) {
        if (where == InKet) {
          expr_ += Scalar(b[dir]) * (Scalar("oo2z") * ABm1 -
                                     Scalar("rho12_over_alpha2") * S_ABm1);
          nflops_ += 5;
        } else {
          expr_ += Scalar(b[dir]) * Scalar("oo2z") * ABm1;
          nflops_ += 3;
        }
      }
    }

    {
      auto S_AB_target = where == InBra ? overlap_factory.make_child(a + _1, b)
                                        : overlap_factory.make_child(a, b + _1);
      if (is_simple()) {
        expr_ += Scalar("two_rho12") * S_AB_target;
        nflops_ += 2;
      }
    }
  }

  // if got here, can decrement by at least 1 quantum
  // add additional derivative terms
  if (deriv) {
    assert(false);  // not yet implemented

    // bf quantum on the build center subtracted by 1
    F a(where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    F b(where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));

    // treatment of derivative terms differs for shell sets and integrals
    // since in computing shell sets transfer/build will occur in all 3
    // directions change in up to all three derivative indices will occur
    for (unsigned int dxyz = 0; dxyz < 3; ++dxyz) {
      if (is_simple() && dxyz != dir)  // for integrals only consider
                                       // derivatives in THE build direction
        continue;

      OriginDerivative<3u> _d1;
      _d1.inc(dxyz);

      std::shared_ptr<DGVertex> _nullptr;

      // dA - _1?
      {
        const OriginDerivative<3u> dAm1(dA - _d1);
        if (exists(dAm1)) {  // yes
          a.deriv() = dAm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ -= Vector(dA)[dxyz] * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ += Vector(dA)[dxyz] * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          a.deriv() = dA;
        }
      }

      // dB - _1?
      {
        const OriginDerivative<3u> dBm1(dB - _d1);
        if (exists(dBm1)) {  // yes
          b.deriv() = dBm1;
          auto AB = factory.make_child(a, b);
          if (is_simple()) {
            if (where == InBra) {  // building on A
              expr_ += Vector(dB)[dxyz] * Scalar("rho12_over_alpha1") * AB;
              nflops_ += 3;
            }
            if (where == InKet) {  // building on B
              expr_ -= Vector(dB)[dxyz] * Scalar("rho12_over_alpha2") * AB;
              nflops_ += 3;
            }
          }
          b.deriv() = dB;
        }
      }
    }
  }  // end of deriv

  return;
}

/** VRR Recurrence Relation for 1-e electrostatic potential integrals.
 * \tparam where specifies whether angular momentum is decreased, in bra or ket.
 */
template <class BFSet, FunctionPosition where>
class VRR_1_ElecPot_1 : public GenericRecurrenceRelation<
                            VRR_1_ElecPot_1<BFSet, where>, BFSet,
                            GenIntegralSet_1_1<BFSet, ElecPotOper, mType> > {
 public:
  typedef VRR_1_ElecPot_1 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_1_1<BFSet, ElecPotOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 9;

  using ParentType::Instance;

  /// Default directionality
  static bool directional() { return ParentType::default_directional(); }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that mainains
  /// registry of these objects
  VRR_1_ElecPot_1(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "OSVRRElecPot"; }
  /** Re-Implementation of GenericRecurrenceRelation::generate_label():
      ElecPot VRR recurrence relation is invariant of m, hence
      to avoid generating identical code make sure that the (unique) label has
     m=0. */
  std::string generate_label() const override {
    typedef typename TargetType::AuxIndexType mType;
    static std::shared_ptr<mType> aux0(new mType(0u));
    std::ostringstream os;
    os << descr() << to_string(where)
       << genintegralset_label(target_->bra(), target_->ket(), aux0,
                               target_->oper());
    return os.str();
  }
};

template <class F, FunctionPosition where>
VRR_1_ElecPot_1<F, where>::VRR_1_ElecPot_1(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  const unsigned int m = Tint->aux()->elem(0);
  const F& _1 = unit<F>(dir);

  {  // can't apply to contracted basis functions
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    if (a.contracted() || b.contracted()) return;
  }

  // if derivative integrals, there will be extra terms (Eq. (143) in Obara &
  // Saika JCP 89 1540 (1988))
  const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero();

  typedef TargetType ChildType;
  ChildFactory<ThisType, ChildType> factory(this);

  // Build on A or B
  {
    // bf quantum on the build center subtracted by 1
    auto a = (where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    if (!exists(a)) return;
    auto b = (where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
    if (!exists(b)) return;

    auto AB_m = factory.make_child(a, b, m);
    auto AB_mp1 = factory.make_child(a, b, m + 1);
    if (is_simple()) {
      expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB_m -
              Vector("PC")[dir] * AB_mp1;
      nflops_ += 3;
    }

    auto am1 = a - _1;
    if (exists(am1)) {
      auto Am1B_m = factory.make_child(am1, b, m);
      auto Am1B_mp1 = factory.make_child(am1, b, m + 1);
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2z") * (Am1B_m - Am1B_mp1);
        nflops_ += 4;
      }
    }
    auto bm1 = b - _1;
    if (exists(bm1)) {
      auto ABm1_m = factory.make_child(a, bm1, m);
      auto ABm1_mp1 = factory.make_child(a, bm1, m + 1);
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2z") * (ABm1_m - ABm1_mp1);
        nflops_ += 4;
      }
    }
  }

  // if got here, can decrement by at least 1 quantum
  // add additional derivative terms
  if (deriv) {
    // bf quantum on the build center subtracted by 1
    F a(where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    F b(where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));

    // treatment of derivative terms differs for shell sets and integrals
    // since in computing shell sets transfer/build will occur in all 3
    // directions change in up to all three derivative indices will occur
    for (unsigned int dxyz = 0; dxyz < 3; ++dxyz) {
      if (is_simple() && dxyz != dir)  // for integrals only consider
                                       // derivatives in THE build direction
        continue;

      OriginDerivative<3u> _d1;
      _d1.inc(dxyz);

      std::shared_ptr<DGVertex> _nullptr;

      // dA - _1?
      {
        const OriginDerivative<3u> dAm1(dA - _d1);
        if (exists(dAm1)) {  // yes
          a.deriv() = dAm1;
          auto AB_m = factory.make_child(a, b, m);
          auto AB_mp1 = factory.make_child(a, b, m + 1);
          if (is_simple()) {
            if (where == InBra) {  // building on A -> derivative of (PA)_i and
                                   // (PC)_i prefactors w.r.t A_i
              expr_ -=
                  Vector(dA)[dxyz] * (Scalar("rho12_over_alpha1") * AB_m +
                                      Scalar("rho12_over_alpha2") * AB_mp1);
              nflops_ += 5;
            }
            if (where == InKet) {  // building on B -> derivative of (PB)_i and
                                   // (PC)_i prefactors w.r.t A_i
              expr_ += Vector(dA)[dxyz] * Scalar("rho12_over_alpha2") *
                       (AB_m - AB_mp1);
              nflops_ += 4;
            }
          }
          a.deriv() = dA;
        }
      }

      // dB - _1?
      {
        const OriginDerivative<3u> dBm1(dB - _d1);
        if (exists(dBm1)) {  // yes
          b.deriv() = dBm1;
          auto AB_m = factory.make_child(a, b, m);
          auto AB_mp1 = factory.make_child(a, b, m + 1);
          if (is_simple()) {
            if (where == InBra) {  // building on A -> derivative of (PA)_i and
                                   // (PC)_i prefactors w.r.t B_i
              expr_ += Vector(dB)[dxyz] * Scalar("rho12_over_alpha1") *
                       (AB_m - AB_mp1);
              nflops_ += 4;
            }
            if (where == InKet) {  // building on B -> derivative of (PB)_i and
                                   // (PC)_i prefactors w.r.t B_i
              expr_ -=
                  Vector(dB)[dxyz] * (Scalar("rho12_over_alpha2") * AB_m +
                                      Scalar("rho12_over_alpha1") * AB_mp1);
              nflops_ += 5;
            }
          }
          b.deriv() = dB;
        }
      }
    }
  }  // end of deriv

  return;
}

/** VRR Recurrence Relation for 1-e spherical multipole moment aka regular solid
 * harmonics integrals. \tparam where specifies whether angular momentum is
 * decreased, in bra or ket.
 */
template <class BFSet, FunctionPosition where>
class VRR_1_SMultipole_1
    : public GenericRecurrenceRelation<
          VRR_1_SMultipole_1<BFSet, where>, BFSet,
          GenIntegralSet_1_1<BFSet, SphericalMultipoleOper, EmptySet> > {
 public:
  typedef VRR_1_SMultipole_1 ThisType;
  typedef BFSet BasisFunctionType;
  typedef SphericalMultipoleOper OperType;
  typedef GenIntegralSet_1_1<BFSet, SphericalMultipoleOper, EmptySet>
      TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 8;

  using ParentType::Instance;

  /// Default directionality
  static bool directional() { return ParentType::default_directional(); }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;
  using typename ParentType::RecurrenceRelation::ExprType;

  static std::string descr() { return "OSVRRSMultipole"; }

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  VRR_1_SMultipole_1(const std::shared_ptr<TargetType>& Tint, unsigned int dir)
      : ParentType(Tint, dir) {
    using namespace libint2::algebra;
    using namespace libint2::prefactor;
    using namespace libint2::braket;
    using Sign = SphericalMultipoleQuanta::Sign;
    using F = BFSet;
    const F& _1 = unit<F>(dir);

    {  // can't apply to contracted basis functions
      F a(Tint->bra(0, 0));
      F b(Tint->ket(0, 0));
      if (a.contracted() || b.contracted()) return;
    }

    // implementation follows J.M. Pérez-Jordá and W. Yang, J Chem Phys 107,
    // 1218 (1997). Eqs. (58-61) for gaussian quanta, and J.M. Pérez-Jordá and
    // W. Yang, J Chem Phys 104, 8003 (1996), Eqs. (23-27) for multipole quanta
    const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
    const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
    const bool deriv = !dA.zero() || !dB.zero();
    if (deriv)  // derivative relations not yet implemented
      return;

    auto O_l_m = Tint->oper()->descr();
    const auto l = O_l_m.l();
    const auto m = O_l_m.m();
    const auto sign = O_l_m.sign();

    typedef TargetType ChildType;
    ChildFactory<ThisType, ChildType> factory(this);

    auto make_child =
        [&](const F& bra, const F& ket,
            const OperType::Descriptor& descr) -> std::shared_ptr<DGVertex> {
      return factory.make_child(bra, ket, EmptySet(), descr);
    };

    // Build on A or B if either bra or ket has quanta, otherwise build
    // multipole quanta
    if (Tint->bra(0, 0).norm() != 0 || Tint->ket(0, 0).norm() != 0) {
      // build A or B

      // bf quantum on the build center subtracted by 1
      auto a = (where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
      if (!exists(a)) return;
      auto b = (where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
      if (!exists(b)) return;

      //
      // first 2 terms are common to Eqs. (58-60)
      //
      auto AB = make_child(a, b, O_l_m);
      if (is_simple()) {
        expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB;
        nflops_ += 1;
      }

      auto am1 = a - _1;
      auto am1_exists = exists(am1);
      auto bm1 = b - _1;
      auto bm1_exists = exists(bm1);

      std::shared_ptr<ExprType> subexpr;  // to be multiplied by 1/(2 zeta)

      if (am1_exists) {
        auto Am1B = make_child(am1, b, O_l_m);
        if (is_simple()) {
          subexpr += Scalar(a[dir]) * Am1B;
          nflops_ += 1;
        }
      }
      if (bm1_exists) {
        auto ABm1 = make_child(a, bm1, O_l_m);
        if (is_simple()) {
          subexpr += Scalar(b[dir]) * ABm1;
          nflops_ += 1;
        }
      }

      // Eq. (58)
      // (a|N_{l,m}|b) += 0.5 * ((a|N_{l-1,m+1}|b) - (a|N_{l-1,m-1}|b))
      auto O_lm1_mp1 = SphericalMultipole_Descr(l - 1, m + 1, sign);
      if (exists(O_lm1_mp1)) {
        auto AB_O_lm1_mp1 = make_child(a, b, O_lm1_mp1);
        if (is_simple() && dir == 0) {
          subexpr += Scalar(O_lm1_mp1.phase() > 0 ? 0.5 : -0.5) * AB_O_lm1_mp1;
          nflops_ += 1;
        }
      }

      auto O_lm1_mm1 = SphericalMultipole_Descr(l - 1, m - 1, sign);
      if (exists(O_lm1_mm1)) {
        auto AB_O_lm1_mm1 = make_child(a, b, O_lm1_mm1);
        if (is_simple() && dir == 0) {
          subexpr -= Scalar(O_lm1_mm1.phase() > 0 ? 0.5 : -0.5) * AB_O_lm1_mm1;
          nflops_ += 1;
        }
      }

      // Eq. (59)
      // (a|N^{+-}_{l,m}|b) {+-}= 0.5 * ((a|N^{-+}_{l-1,m+1}|b) +
      // (a|N^{-+}_{l-1,m-1}|b))
      const auto m_pfac = sign == Sign::plus ? 1 : -1;
      auto Om_lm1_mp1 = SphericalMultipole_Descr(l - 1, m + 1, flip(sign));
      if (exists(Om_lm1_mp1)) {
        auto AB_Om_lm1_mp1 = make_child(a, b, Om_lm1_mp1);
        if (is_simple() && dir == 1) {
          subexpr += Scalar(m_pfac * (Om_lm1_mp1.phase() > 0 ? 0.5 : -0.5)) *
                     AB_Om_lm1_mp1;
          nflops_ += 1;
        }
      }

      auto Om_lm1_mm1 = SphericalMultipole_Descr(l - 1, m - 1, flip(sign));
      if (exists(Om_lm1_mm1)) {
        auto AB_Om_lm1_mm1 = make_child(a, b, Om_lm1_mm1);
        if (is_simple() && dir == 1) {
          subexpr += Scalar(m_pfac * (Om_lm1_mm1.phase() > 0 ? 0.5 : -0.5)) *
                     AB_Om_lm1_mm1;
          nflops_ += 1;
        }
      }

      // Eq. (60)
      auto O_lm1_m = SphericalMultipole_Descr(l - 1, m, sign);
      if (exists(O_lm1_m)) {
        auto AB_O_lm1_m = make_child(a, b, O_lm1_m);
        if (is_simple() && dir == 2) {
          subexpr += AB_O_lm1_m;
          nflops_ += 1;
        }
      }

      if (is_simple()) {
        if (subexpr) expr_ += Scalar("oo2z") * subexpr;
      }
    } else {
      // build multipole quanta only if dir == 0
      if (dir != 0) return;
      auto a = Tint->bra(0, 0);
      auto b = Tint->ket(0, 0);

      if (l == m) {
        if (l == 0) {                  // Eq.
          assert(sign == Sign::plus);  // Eq. (24) should not be needed
                                       // since N^-_{0,0} should just be
                                       // omitted above
          typedef GenIntegralSet_1_1<BFSet, CartesianMultipoleOper<3u>,
                                     EmptySet>
              OverlapType;
          ChildFactory<ThisType, OverlapType> overlap_factory(this);
          auto S00 = overlap_factory.make_child(a, b);
          if (is_simple()) {
            expr_ = Scalar(1) * S00;  // Eq. (25) and Eq. (61)
          }
        } else {
          std::shared_ptr<ExprType> subexpr;  // to be multiplied by - 1/(2 m)

          // Eqs. (25-26)
          // (0|N^+_{m,m}|0) = -(1/(2m)) ( x (0|N^+_{m-1,m-1}|0) - y
          // (0|N^-_{m-1,m-1}|0) )
          auto Op_lm1_mm1 = SphericalMultipole_Descr(l - 1, m - 1, Sign::plus);
          auto AB_Op_lm1_mm1 = make_child(a, b, Op_lm1_mm1);
          if (is_simple()) {
            subexpr =
                Scalar(sign == Sign::plus ? "PO_x" : "PO_y") * AB_Op_lm1_mm1;
          }
          auto Om_lm1_mm1 = SphericalMultipole_Descr(l - 1, m - 1, Sign::minus);
          if (exists(Om_lm1_mm1)) {
            auto AB_Om_lm1_mm1 = make_child(a, b, Om_lm1_mm1);
            if (is_simple()) {
              if (sign == Sign::plus)
                subexpr -= Scalar("PO_y") * AB_Om_lm1_mm1;
              else  // sign == Sign::minus
                subexpr += Scalar("PO_x") * AB_Om_lm1_mm1;
            }
          }
          if (is_simple()) {
            assert(subexpr);
            expr_ = Scalar(-0.5 / m) * subexpr;
          }
        }
      } else {  // l != m, use Eq. 27
        std::shared_ptr<ExprType>
            subexpr;  // to be multiplied by 1/((l+m)(l-m))

        auto O_lm1_m = SphericalMultipole_Descr(l - 1, m, sign);
        if (exists(O_lm1_m)) {
          auto AB_O_lm1_m = make_child(a, b, O_lm1_m);
          if (is_simple())
            subexpr = Scalar((2 * l - 1)) * Scalar("PO_z") * AB_O_lm1_m;
        }

        auto O_lm2_m = SphericalMultipole_Descr(l - 2, m, sign);
        if (exists(O_lm2_m)) {
          auto AB_O_lm2_m = make_child(a, b, O_lm2_m);
          if (is_simple()) subexpr -= Scalar("PO2") * AB_O_lm2_m;
        }

        if (is_simple()) {
          assert(subexpr);
          expr_ = Scalar(1.0 / ((l + m) * (l - m))) * subexpr;
        }
      }
    }

    return;
  }

};  // VRR_1_SMultipole_1

};  // namespace libint2

#endif
