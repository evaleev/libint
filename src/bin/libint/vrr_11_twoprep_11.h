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

#ifndef _libint2_src_bin_libint_vrr11twoprep11_h_
#define _libint2_src_bin_libint_vrr11twoprep11_h_

#include <generic_rr.h>
#include <twoprep_11_11.h>

namespace libint2 {

/** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
the angular momentum is raised. where specifies whether angular momentum is
decreased in bra or ket.
*/
template <class BFSet, int part, FunctionPosition where>
class VRR_11_TwoPRep_11 : public GenericRecurrenceRelation<
                              VRR_11_TwoPRep_11<BFSet, part, where>, BFSet,
                              GenIntegralSet_11_11<BFSet, TwoPRep, mType> > {
 public:
  typedef VRR_11_TwoPRep_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, TwoPRep, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 26;

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
  VRR_11_TwoPRep_11(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "OSVRR"; }
  /** Re-Implementation of GenericRecurrenceRelation::generate_label():
      TwoPRep VRR recurrence relations codes are independent of m (it never
     appears anywhere in equations), hence to avoid generating identical code
     make sure that the (unique) label has m=0. */
  std::string generate_label() const override {
    typedef typename TargetType::AuxIndexType mType;
    static std::shared_ptr<mType> aux0(new mType(0u));
    std::ostringstream os;
    os << descr() << "P" << part << to_string(where)
       << genintegralset_label(target_->bra(), target_->ket(), aux0,
                               target_->oper());
    return os.str();
  }

#if LIBINT_ENABLE_GENERIC_CODE
  /// Implementation of RecurrenceRelation::has_generic()
  bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const override;
  /// Implementation of RecurrenceRelation::generic_header()
  std::string generic_header() const override;
  /// Implementation of RecurrenceRelation::generic_instance()
  std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const override;
#endif
};

template <class F, int part, FunctionPosition where>
VRR_11_TwoPRep_11<F, part, where>::VRR_11_TwoPRep_11(
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
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));
    if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
      return;
  }

  // if derivative integrals, there will be extra terms (Eq. (143) in Obara &
  // Saika JCP 89)
  const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = Tint->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = Tint->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();

  // This is a hack to avoid creating recurrence relations for which generic
  // code has not been yet implemented
#if LIBINT_ENABLE_GENERIC_CODE
  {
    F sh_a(target_->bra(0, 0));
    F sh_b(target_->ket(0, 0));
    F sh_c(target_->bra(1, 0));
    F sh_d(target_->ket(1, 0));
    // generic code works for a0c0 of 0a0c classes where am(a) > 1 and am(c) > 1
    // to generate optimized code for xxxx integral need to generate specialized
    // code for up to (x+x)0(x+x)0 integrals
    if (sh_b.zero() && sh_d.zero() &&
        (sh_a.norm() > 1u &&
         sh_c.norm() > 1u)) {  // have a generic implemented ...
      if (part != 0)  // ... but only implemented build on A in this case
        return;
    }
    if (sh_a.zero() && sh_c.zero() && (sh_b.norm() > 1u && sh_d.norm() > 1u)) {
      if (part != 0)  // but only implemented build on B in this case
        return;
    }
  }
#endif

  typedef TargetType ChildType;
  ChildFactory<ThisType, ChildType> factory(this);

  bool part0_has_unit = false, part1_has_unit = false;

  // Build on A
  if (part == 0 && where == InBra) {
    F a(Tint->bra(0, 0) - _1);
    if (!exists(a)) return;
    F b(Tint->ket(0, 0));
    const bool unit_b = (b == F::unit());
    part0_has_unit |= unit_b;
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));

    std::shared_ptr<DGVertex> ABCD_m;
    if (not unit_b) ABCD_m = factory.make_child(a, b, c, d, m);
    auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
    if (is_simple()) {
      if (not unit_b) {
        expr_ = Vector("PA")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 3;
      } else {
        expr_ = Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 1;
      }
    }

    // simplified 3-center VRR due to Ahlrichs (PCCP 6, 5119 (2004))
    const bool ahlrichs_simplification = a.pure_sh() && unit_b;
    auto am1 = a - _1;
    if (exists(am1) && not ahlrichs_simplification) {
      auto Am1BCD_m = factory.make_child(am1, b, c, d, m);
      auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(a[dir]) * Scalar("oo2z") *
                 (Scalar("roz") * Am1BCD_mp1 - Am1BCD_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2z") *
                 (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& bm1 = b - _1;
    if (exists(bm1)) {
      auto ABm1CD_m = factory.make_child(a, bm1, c, d, m);
      auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(b[dir]) * Scalar("oo2z") *
                 (Scalar("roz") * ABm1CD_mp1 - ABm1CD_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2z") *
                 (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& cm1 = c - _1;
    if (exists(cm1)) {
      auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2ze") * ABCm1D_mp1;
        nflops_ += 3;
      }
    }
    const F& dm1 = d - _1;
    if (exists(dm1)) {
      auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1);
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2ze") * ABCDm1_mp1;
        nflops_ += 3;
      }
    }
  }
  // Build on B
  if (part == 0 && where == InKet) {
    F a(Tint->bra(0, 0));
    const bool unit_a = (a == F::unit());
    part0_has_unit |= unit_a;
    F b(Tint->ket(0, 0) - _1);
    if (!exists(b)) return;
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));

    std::shared_ptr<DGVertex> ABCD_m;
    if (not unit_a) ABCD_m = factory.make_child(a, b, c, d, m);
    auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
    if (is_simple()) {
      if (not unit_a) {
        expr_ = Vector("PB")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 3;
      } else {
        expr_ = Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 1;
      }
    }

    const F& am1 = a - _1;
    if (exists(am1)) {
      auto Am1BCD_m = factory.make_child(am1, b, c, d, m);
      auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(a[dir]) * Scalar("oo2z") *
                 (Scalar("roz") * Am1BCD_mp1 - Am1BCD_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2z") *
                 (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);
        nflops_ += 5;
      }
#endif
    }
    // simplified 3-center VRR due to Ahlrichs (PCCP 6, 5119 (2004))
    const bool ahlrichs_simplification = b.pure_sh() && unit_a;
    const F& bm1 = b - _1;
    if (exists(bm1) && not ahlrichs_simplification) {
      auto ABm1CD_m = factory.make_child(a, bm1, c, d, m);
      auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(b[dir]) * Scalar("oo2z") *
                 (Scalar("roz") * ABm1CD_mp1 - ABm1CD_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2z") *
                 (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& cm1 = c - _1;
    if (exists(cm1)) {
      auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2ze") * ABCm1D_mp1;
        nflops_ += 3;
      }
    }
    const F& dm1 = d - _1;
    if (exists(dm1)) {
      auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1);
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2ze") * ABCDm1_mp1;
        nflops_ += 3;
      }
    }
  }
  // Build on C
  if (part == 1 && where == InBra) {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0) - _1);
    if (!exists(c)) return;
    F d(Tint->ket(1, 0));
    const bool unit_d = (d == F::unit());
    part1_has_unit |= unit_d;

    std::shared_ptr<DGVertex> ABCD_m;
    if (not unit_d) ABCD_m = factory.make_child(a, b, c, d, m);
    auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
    if (is_simple()) {
      if (not unit_d) {
        expr_ = Vector("QC")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 3;
      } else {
        expr_ = Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 1;
      }
    }

    // simplified 3-center VRR due to Ahlrichs (PCCP 6, 5119 (2004))
    const bool ahlrichs_simplification = c.pure_sh() && unit_d;
    const F& cm1 = c - _1;
    if (exists(cm1) && not ahlrichs_simplification) {
      auto ABCm1D_m = factory.make_child(a, b, cm1, d, m);
      auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(c[dir]) * Scalar("oo2e") *
                 (Scalar("roe") * ABCm1D_mp1 - ABCm1D_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2e") *
                 (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& dm1 = d - _1;
    if (exists(dm1)) {
      auto ABCDm1_m = factory.make_child(a, b, c, dm1, m);
      auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(d[dir]) * Scalar("oo2e") *
                 (Scalar("roe") * ABCDm1_mp1 - ABCDm1_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2e") *
                 (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& am1 = a - _1;
    if (exists(am1)) {
      auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2ze") * Am1BCD_mp1;
        nflops_ += 3;
      }
    }
    const F& bm1 = b - _1;
    if (exists(bm1)) {
      auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2ze") * ABm1CD_mp1;
        nflops_ += 3;
      }
    }
  }
  // Build on D
  if (part == 1 && where == InKet) {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0));
    const bool unit_c = (c == F::unit());
    part1_has_unit |= unit_c;
    F d(Tint->ket(1, 0) - _1);
    if (!exists(d)) return;

    std::shared_ptr<DGVertex> ABCD_m;
    if (not unit_c) ABCD_m = factory.make_child(a, b, c, d, m);
    auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
    if (is_simple()) {
      if (not unit_c) {
        expr_ = Vector("QD")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 3;
      } else {
        expr_ = Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 1;
      }
    }

    const F& cm1 = c - _1;
    if (exists(cm1)) {
      auto ABCm1D_m = factory.make_child(a, b, cm1, d, m);
      auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(c[dir]) * Scalar("oo2e") *
                 (Scalar("roe") * ABCm1D_mp1 - ABCm1D_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(c[dir]) * Scalar("oo2e") *
                 (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);
        nflops_ += 5;
      }
#endif
    }
    // simplified 3-center VRR due to Ahlrichs (PCCP 6, 5119 (2004))
    const bool ahlrichs_simplification = d.pure_sh() && unit_c;
    const F& dm1 = d - _1;
    if (exists(dm1) && not ahlrichs_simplification) {
      auto ABCDm1_m = factory.make_child(a, b, c, dm1, m);
      auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1);
#if LIBINT_GENERATE_FMA
      // this form is amenable to generation of fmsub
      if (is_simple()) {
        expr_ -= Scalar(d[dir]) * Scalar("oo2e") *
                 (Scalar("roe") * ABCDm1_mp1 - ABCDm1_m);
        nflops_ += 5;
      }
#else
      if (is_simple()) {
        expr_ += Scalar(d[dir]) * Scalar("oo2e") *
                 (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);
        nflops_ += 5;
      }
#endif
    }
    const F& am1 = a - _1;
    if (exists(am1)) {
      auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(a[dir]) * Scalar("oo2ze") * Am1BCD_mp1;
        nflops_ += 3;
      }
    }
    const F& bm1 = b - _1;
    if (exists(bm1)) {
      auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1);
      if (is_simple()) {
        expr_ += Scalar(b[dir]) * Scalar("oo2ze") * ABm1CD_mp1;
        nflops_ += 3;
      }
    }
  }

  // if got here, can decrement by at least 1 quantum
  // add additional derivative terms
  if (deriv) {
    // bf quantum on the build center subtracted by 1
    F a(part == 0 && where == InBra ? Tint->bra(0, 0) - _1 : Tint->bra(0, 0));
    F b(part == 0 && where == InKet ? Tint->ket(0, 0) - _1 : Tint->ket(0, 0));
    F c(part == 1 && where == InBra ? Tint->bra(1, 0) - _1 : Tint->bra(1, 0));
    F d(part == 1 && where == InKet ? Tint->ket(1, 0) - _1 : Tint->ket(1, 0));

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
          auto ABCD_m = (part == 0 && not part0_has_unit)
                            ? factory.make_child(a, b, c, d, m)
                            : _nullptr;
          auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
          if (is_simple()) {
            if (part == 0 && where == InBra) {  // building on A
              if (not part0_has_unit) {
                expr_ -= Vector(dA)[dxyz] *
                         (Scalar("rho12_over_alpha1") * ABCD_m +
                          Scalar("alpha1_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dA)[dxyz] *
                         (Scalar("alpha1_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 0 && where == InKet) {  // building on B
              if (not part0_has_unit) {
                expr_ += Vector(dA)[dxyz] *
                         (Scalar("rho12_over_alpha2") * ABCD_m -
                          Scalar("alpha1_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dA)[dxyz] *
                         (Scalar("alpha1_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 1) {  // building on C or D
              expr_ += Vector(dA)[dxyz] * Scalar("alpha1_over_zetapluseta") *
                       ABCD_mp1;
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
          auto ABCD_m = (part == 0 && not part0_has_unit)
                            ? factory.make_child(a, b, c, d, m)
                            : _nullptr;
          auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
          if (is_simple()) {
            if (part == 0 && where == InBra) {  // building on A
              if (not part0_has_unit) {
                expr_ += Vector(dB)[dxyz] *
                         (Scalar("rho12_over_alpha1") * ABCD_m -
                          Scalar("alpha2_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dB)[dxyz] *
                         (Scalar("alpha2_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 0 && where == InKet) {  // building on B
              if (not part0_has_unit) {
                expr_ -= Vector(dB)[dxyz] *
                         (Scalar("rho12_over_alpha2") * ABCD_m +
                          Scalar("alpha2_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dB)[dxyz] *
                         (Scalar("alpha2_rho_over_zeta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 1) {  // building on C or D
              expr_ += Vector(dB)[dxyz] * Scalar("alpha2_over_zetapluseta") *
                       ABCD_mp1;
              nflops_ += 3;
            }
          }
          b.deriv() = dB;
        }
      }

      // dC - _1?
      {
        const OriginDerivative<3u> dCm1(dC - _d1);
        if (exists(dCm1)) {  // yes
          c.deriv() = dCm1;
          auto ABCD_m = (part == 1 && not part1_has_unit)
                            ? factory.make_child(a, b, c, d, m)
                            : _nullptr;
          auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
          if (is_simple()) {
            if (part == 0) {  // building on A or B
              expr_ += Vector(dC)[dxyz] * Scalar("alpha3_over_zetapluseta") *
                       ABCD_mp1;
              nflops_ += 3;
            }
            if (part == 1 && where == InBra) {  // building on C
              if (not part1_has_unit) {
                expr_ -= Vector(dC)[dxyz] *
                         (Scalar("rho34_over_alpha3") * ABCD_m +
                          Scalar("alpha3_rho_over_eta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dC)[dxyz] *
                         (Scalar("alpha3_rho_over_eta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 1 && where == InKet) {  // building on D
              if (not part1_has_unit) {
                expr_ += Vector(dC)[dxyz] *
                         (Scalar("rho34_over_alpha4") * ABCD_m -
                          Scalar("alpha3_rho_over_eta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dC)[dxyz] *
                         (Scalar("alpha3_rho_over_eta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
          }
          c.deriv() = dC;
        }
      }

      // dD - _1?
      {
        const OriginDerivative<3u> dDm1(dD - _d1);
        if (exists(dDm1)) {  // yes
          d.deriv() = dDm1;
          auto ABCD_m = (part == 1 && not part1_has_unit)
                            ? factory.make_child(a, b, c, d, m)
                            : _nullptr;
          auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1);
          if (is_simple()) {
            if (part == 0) {  // building on A or B
              expr_ += Vector(dD)[dxyz] * Scalar("alpha4_over_zetapluseta") *
                       ABCD_mp1;
              nflops_ += 3;
            }
            if (part == 1 && where == InBra) {  // building on C
              if (not part1_has_unit) {
                expr_ += Vector(dD)[dxyz] *
                         (Scalar("rho34_over_alpha3") * ABCD_m -
                          Scalar("alpha4_rho_over_eta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dD)[dxyz] *
                         (Scalar("alpha4_rho_over_eta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
            if (part == 1 && where == InKet) {  // building on D
              if (not part1_has_unit) {
                expr_ -= Vector(dD)[dxyz] *
                         (Scalar("rho34_over_alpha4") * ABCD_m +
                          Scalar("alpha4_rho_over_eta2") * ABCD_mp1);
                nflops_ += 5;
              } else {
                expr_ -= Vector(dD)[dxyz] *
                         (Scalar("alpha4_rho_over_eta2") * ABCD_mp1);
                nflops_ += 3;
              }
            }
          }
          d.deriv() = dD;
        }
      }
    }
  }  // end of deriv

  return;
}

#if LIBINT_ENABLE_GENERIC_CODE
template <class F, int part, FunctionPosition where>
bool VRR_11_TwoPRep_11<F, part, where>::has_generic(
    const std::shared_ptr<CompilationParameters>& cparams) const {
  if (TrivialBFSet<F>::result) return false;

  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const unsigned int max_opt_am = cparams->max_am_opt();
  // generic code works for a0c0 of 0a0c classes where am(a) > 1 and am(c) > 1
  // to generate optimized code for xxxx integral need to generate specialized
  // code for up to (x+x)0(x+x)0 integrals
  if (sh_b.zero() && sh_d.zero() &&
      (sh_a.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_c.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_a.norm() > 1u && sh_c.norm() > 1u)) {
    assert(part == 0);  // has only implemented build on A in this case
    return true;
  }
  if (sh_a.zero() && sh_c.zero() &&
      (sh_b.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_d.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_b.norm() > 1u && sh_d.norm() > 1u)) {
    assert(part == 0);  // has only implemented build on B in this case
    return true;
  }
  return false;
}

template <class F, int part, FunctionPosition where>
std::string VRR_11_TwoPRep_11<F, part, where>::generic_header() const {
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();

  const OriginDerivative<3u> dA = target_->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = target_->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = target_->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = target_->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();

  if (!deriv) {
    if (xsxs) return std::string("OSVRR_xs_xs.h");
    if (sxsx) return std::string("OSVRR_sx_sx.h");
  } else {
    if (xsxs) return std::string("OSVRR_xs_xs_deriv.h");
    if (sxsx) return std::string("OSVRR_sx_sx_deriv.h");
  }
  abort();  // unreachable
}

template <class F, int part, FunctionPosition where>
std::string VRR_11_TwoPRep_11<F, part, where>::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  std::ostringstream oss;
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();
  bool ahlrichs_simplification = false;
  bool unit_s = false;
  bool part0_has_unit = false;
  bool part1_has_unit = false;
  if (xsxs) {
    ahlrichs_simplification = (sh_a.pure_sh() && sh_b.is_unit()) ||
                              (sh_c.pure_sh() && sh_d.is_unit());
    unit_s = sh_b.is_unit();
    part0_has_unit = sh_b.is_unit();
    part1_has_unit = sh_d.is_unit();
  }
  if (sxsx) {
    ahlrichs_simplification = (sh_b.pure_sh() && sh_a.is_unit()) ||
                              (sh_d.pure_sh() && sh_c.is_unit());
    unit_s = sh_a.is_unit();
    part0_has_unit = sh_a.is_unit();
    part1_has_unit = sh_c.is_unit();
  }

  const OriginDerivative<3u> dA = target_->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = target_->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = target_->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = target_->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();

  oss << "using namespace libint2;" << std::endl;

  if (!deriv) {  // for regular integrals I know exactly how many prerequisites
                 // I need
    if (xsxs) {
      oss << "libint2::OS" << (ahlrichs_simplification ? "A" : "")
          << "VRR_xs_xs<" << part << "," << sh_a.norm() << "," << sh_c.norm()
          << ",";
    }
    if (sxsx) {
      oss << "libint2::OS" << (ahlrichs_simplification ? "A" : "")
          << "VRR_sx_sx<" << part << "," << sh_b.norm() << "," << sh_d.norm()
          << ",";
    }
    if (not ahlrichs_simplification) oss << (unit_s ? "true," : "false,");
    oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
    oss << ">::compute(inteval";

    oss << "," << args->symbol(0);  // target
    if (not ahlrichs_simplification &&
        unit_s)     // purely to avoid having a 4-term generic RR, reuse 5-term
                    // with dummy argument
      oss << ",0";  // src0-> 0x0
    const unsigned int nargs = args->n();
    for (unsigned int a = 1; a < nargs; a++) {
      oss << "," << args->symbol(a);
    }
    oss << ");";
  } else {  // deriv == true -> only some arguments are needed
    if (xsxs) {
      oss << "libint2::OS" << (ahlrichs_simplification ? "A" : "")
          << "VRR_xs_xs_deriv<" << part << "," << sh_a.norm() << ","
          << sh_c.norm() << ",";
    }
    if (sxsx) {
      oss << "libint2::OS" << (ahlrichs_simplification ? "A" : "")
          << "VRR_sx_sx_deriv<" << part << "," << sh_b.norm() << ","
          << sh_d.norm() << ",";
    }

    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      oss << sh_a.deriv().d(xyz) << ",";
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      oss << sh_b.deriv().d(xyz) << ",";
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      oss << sh_c.deriv().d(xyz) << ",";
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      oss << sh_d.deriv().d(xyz) << ",";
    if (not ahlrichs_simplification) oss << (unit_s ? "true," : "false,");
    oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
    oss << ">::compute(inteval";
    // out of all 22 possible prerequisites first few have same derivative
    // degree as the target 5 if standard 4-center integral 1+4 if the center
    // opposite the build center carries a unit function
    //     will pass 0 instead of src0
    // 2 if Ahlrichs
    oss << "," << args->symbol(0);  // target
    if (not ahlrichs_simplification &&
        unit_s)     // purely to avoid having a 4-term generic RR, reuse 5-term
                    // with dummy argument
      oss << ",0";  // src0-> 0x0
    // const unsigned int nargs = args->n();
    unsigned int arg = 1;
    for (;
         arg <
         (ahlrichs_simplification ? 3 : (unit_s ? 5 : 6));  // nargs + 1 target
         arg++) {
      oss << "," << args->symbol(arg);
    }
    for (unsigned int xyz = 0; xyz < 3; ++xyz) {
      if (sh_a.deriv().d(xyz) > 0) {
        // see the dA-1 clause: (ab|cd)^m is skipped
        oss << ","
            << (part0_has_unit ? std::to_string(0) : args->symbol(arg++));
        oss << "," << args->symbol(arg++);
      } else
        oss << ",0,0";
      if (sh_b.deriv().d(xyz) > 0) {
        // see the dB-1 clause: (ab|cd)^m is skipped
        oss << ","
            << (part0_has_unit ? std::to_string(0) : args->symbol(arg++));
        oss << "," << args->symbol(arg++);
      } else
        oss << ",0,0";
      if (sh_c.deriv().d(xyz) > 0) {
        oss << "," << args->symbol(arg++);
      } else
        oss << ",0";
      if (sh_d.deriv().d(xyz) > 0) {
        oss << "," << args->symbol(arg++);
      } else
        oss << ",0";
    }
    oss << ");";
  }

  return oss.str();
}
#endif  // #if !LIBINT_ENABLE_GENERIC_CODE

};  // namespace libint2

#endif
