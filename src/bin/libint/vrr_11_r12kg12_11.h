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

#ifndef _libint2_src_bin_libint_vrr11r12kg1211_h_
#define _libint2_src_bin_libint_vrr11r12kg1211_h_

#include <generic_rr.h>
#include <r12kg12_11_11.h>

namespace libint2 {

/** VRR Recurrence Relation for 2-e integrals of the R12_k_G12 operators.
part specifies the angular momentum of which particle is raised.
where specifies whether the angular momentum is decreased in bra or ket.
*/
template <class BFSet, int part, FunctionPosition where>
class VRR_11_R12kG12_11 : public GenericRecurrenceRelation<
                              VRR_11_R12kG12_11<BFSet, part, where>, BFSet,
                              GenIntegralSet_11_11<BFSet, R12kG12, mType> > {
 public:
  typedef VRR_11_R12kG12_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, R12kG12, mType> TargetType;
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

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  VRR_11_R12kG12_11(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "VRR"; }
  /** Re-Implementation of GenericRecurrenceRelation::generate_label():
      R12kG12 VRR recurrence relations codes are independent of m (it never
     appears anywhere in equations), hence to avoid generating identical code
     make sure that the (unique) label does not contain m. */
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
VRR_11_R12kG12_11<F, part, where>::VRR_11_R12kG12_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  const int K = Tint->oper()->descr().K();
  typedef R12_k_G12_Descr R12kG12Descr;
  const R12kG12 oK((R12kG12Descr(K)));
  const unsigned int m = Tint->aux()->elem(0);
  const F _1 = unit<F>(dir);

  // does not work for contracted functions
  {
    F a(Tint->bra(0, 0));
    F b(Tint->ket(0, 0));
    F c(Tint->bra(1, 0));
    F d(Tint->ket(1, 0));
    if (a.contracted() || b.contracted() || c.contracted() || d.contracted() ||
        Tint->oper()->descr().contracted())
      return;
  }

  // does not work for derivative integrals (yet or ever)
  const OriginDerivative<3u> dA = Tint->bra(0, 0).deriv();
  const OriginDerivative<3u> dB = Tint->ket(0, 0).deriv();
  const OriginDerivative<3u> dC = Tint->bra(1, 0).deriv();
  const OriginDerivative<3u> dD = Tint->ket(1, 0).deriv();
  const bool deriv = !dA.zero() || !dB.zero() || !dC.zero() || !dD.zero();
  assert(!deriv);

  typedef TargetType ChildType;
  ChildFactory<ThisType, ChildType> factory(this);

  // if K is -1, the recurrence relation looks exactly as it would for ERI
  // thus generate the same code, and remember to use appropriate prefactors
  if (K == -1) {
    // Build on A
    if (part == 0 && where == InBra) {
      F a(Tint->bra(0, 0) - _1);
      if (!exists(a)) return;
      F b(Tint->ket(0, 0));
      F c(Tint->bra(1, 0));
      F d(Tint->ket(1, 0));

      auto ABCD_m = factory.make_child(a, b, c, d, m, oK);
      auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1, oK);
      if (is_simple()) {
        expr_ = Vector("PA")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 3;
      }

      const F& am1 = a - _1;
      if (exists(am1)) {
        auto Am1BCD_m = factory.make_child(am1, b, c, d, m, oK);
        auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(a[dir]) * Scalar("oo2z") *
                   (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);
          nflops_ += 5;
        }
      }
      const F& bm1 = b - _1;
      if (exists(bm1)) {
        auto ABm1CD_m = factory.make_child(a, bm1, c, d, m, oK);
        auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(b[dir]) * Scalar("oo2z") *
                   (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);
          nflops_ += 5;
        }
      }
      const F& cm1 = c - _1;
      if (exists(cm1)) {
        auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(c[dir]) * Scalar("oo2ze") * ABCm1D_mp1;
          nflops_ += 3;
        }
      }
      const F& dm1 = d - _1;
      if (exists(dm1)) {
        auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(d[dir]) * Scalar("oo2ze") * ABCDm1_mp1;
          nflops_ += 3;
        }
      }
      return;
    }
    // Build on A
    if (part == 0 && where == InKet) {
      F a(Tint->bra(0, 0));
      F b(Tint->ket(0, 0) - _1);
      if (!exists(b)) return;
      F c(Tint->bra(1, 0));
      F d(Tint->ket(1, 0));

      auto ABCD_m = factory.make_child(a, b, c, d, m, oK);
      auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1, oK);
      if (is_simple()) {
        expr_ = Vector("PB")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;
        nflops_ += 3;
      }

      const F& am1 = a - _1;
      if (exists(am1)) {
        auto Am1BCD_m = factory.make_child(am1, b, c, d, m, oK);
        auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(a[dir]) * Scalar("oo2z") *
                   (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);
          nflops_ += 5;
        }
      }
      const F& bm1 = b - _1;
      if (exists(bm1)) {
        auto ABm1CD_m = factory.make_child(a, bm1, c, d, m, oK);
        auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(b[dir]) * Scalar("oo2z") *
                   (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);
          nflops_ += 5;
        }
      }
      const F& cm1 = c - _1;
      if (exists(cm1)) {
        auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(c[dir]) * Scalar("oo2ze") * ABCm1D_mp1;
          nflops_ += 3;
        }
      }
      const F& dm1 = d - _1;
      if (exists(dm1)) {
        auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(d[dir]) * Scalar("oo2ze") * ABCDm1_mp1;
          nflops_ += 3;
        }
      }
      return;
    }
    // Build on C
    if (part == 1 && where == InBra) {
      F a(Tint->bra(0, 0));
      F b(Tint->ket(0, 0));
      F c(Tint->bra(1, 0) - _1);
      if (!exists(c)) return;
      F d(Tint->ket(1, 0));

      auto ABCD_m = factory.make_child(a, b, c, d, m, oK);
      auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1, oK);
      if (is_simple()) {
        expr_ = Vector("QC")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 3;
      }

      const F& cm1 = c - _1;
      if (exists(cm1)) {
        auto ABCm1D_m = factory.make_child(a, b, cm1, d, m, oK);
        auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(c[dir]) * Scalar("oo2e") *
                   (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);
          nflops_ += 5;
        }
      }
      const F& dm1 = d - _1;
      if (exists(dm1)) {
        auto ABCDm1_m = factory.make_child(a, b, c, dm1, m, oK);
        auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(d[dir]) * Scalar("oo2e") *
                   (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);
          nflops_ += 5;
        }
      }
      const F& am1 = a - _1;
      if (exists(am1)) {
        auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(a[dir]) * Scalar("oo2ze") * Am1BCD_mp1;
          nflops_ += 3;
        }
      }
      const F& bm1 = b - _1;
      if (exists(bm1)) {
        auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(b[dir]) * Scalar("oo2ze") * ABm1CD_mp1;
          nflops_ += 3;
        }
      }
      return;
    }
    // Build on D
    if (part == 1 && where == InKet) {
      F a(Tint->bra(0, 0));
      F b(Tint->ket(0, 0));
      F c(Tint->bra(1, 0));
      F d(Tint->ket(1, 0) - _1);
      if (!exists(d)) return;

      auto ABCD_m = factory.make_child(a, b, c, d, m, oK);
      auto ABCD_mp1 = factory.make_child(a, b, c, d, m + 1, oK);
      if (is_simple()) {
        expr_ = Vector("QD")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;
        nflops_ += 3;
      }

      const F& cm1 = c - _1;
      if (exists(cm1)) {
        auto ABCm1D_m = factory.make_child(a, b, cm1, d, m, oK);
        auto ABCm1D_mp1 = factory.make_child(a, b, cm1, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(c[dir]) * Scalar("oo2e") *
                   (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);
          nflops_ += 5;
        }
      }
      const F& dm1 = d - _1;
      if (exists(dm1)) {
        auto ABCDm1_m = factory.make_child(a, b, c, dm1, m, oK);
        auto ABCDm1_mp1 = factory.make_child(a, b, c, dm1, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(d[dir]) * Scalar("oo2e") *
                   (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);
          nflops_ += 5;
        }
      }
      const F& am1 = a - _1;
      if (exists(am1)) {
        auto Am1BCD_mp1 = factory.make_child(am1, b, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(a[dir]) * Scalar("oo2ze") * Am1BCD_mp1;
          nflops_ += 3;
        }
      }
      const F& bm1 = b - _1;
      if (exists(bm1)) {
        auto ABm1CD_mp1 = factory.make_child(a, bm1, c, d, m + 1, oK);
        if (is_simple()) {
          expr_ += Scalar(b[dir]) * Scalar("oo2ze") * ABm1CD_mp1;
          nflops_ += 3;
        }
      }
      return;
    }
    return;
  }  // K == -1?
  else {
    // K != -1, the auxiliary quantum number is not used
    if (m != 0)
      throw std::logic_error(
          "VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- "
          "nonzero auxiliary quantum detected.");

    // can build (a0|c0) or (0b|0d)
    bool xsxs = false;
    bool sxsx = false;
    {
      F b(Tint->ket(0, 0) - _1);
      F d(Tint->ket(1, 0) - _1);
      if (!exists(b) && !exists(d)) xsxs = true;
    }
    {
      F a(Tint->bra(0, 0) - _1);
      F c(Tint->bra(1, 0) - _1);
      if (!exists(a) && !exists(c)) sxsx = true;
    }
    // can't handle the general case
    if (!xsxs && !sxsx) return;
    // can't handle (ss|ss) case
    if (xsxs && sxsx) return;

    if (xsxs) {
      // Build on A
      if (part == 0 && where == InBra) {
        F a(Tint->bra(0, 0) - _1);
        if (!exists(a)) return;
        F b(Tint->ket(0, 0));
        F c(Tint->bra(1, 0));
        F d(Tint->ket(1, 0));

        auto ABCD_K = factory.make_child(a, b, c, d, 0u, oK);
        if (is_simple()) {
          expr_ = Vector("R12kG12_pfac0_0")[dir] * ABCD_K;
          nflops_ += 1;
        }
        const F& am1 = a - _1;
        if (exists(am1)) {
          auto Am1BCD_K = factory.make_child(am1, b, c, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(a[dir]) * Scalar("R12kG12_pfac1_0") * Am1BCD_K;
            nflops_ += 3;
          }
        }
        const F& cm1 = c - _1;
        if (exists(cm1)) {
          auto ABCm1D_K = factory.make_child(a, b, cm1, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(c[dir]) * Scalar("R12kG12_pfac2") * ABCm1D_K;
            nflops_ += 3;
          }
        }
        if (K != 0) {
          const R12kG12 oKm2((R12kG12Descr(K - 2)));
          auto Ap1BCD_Km2 = factory.make_child(a + _1, b, c, d, 0u, oKm2);
          auto ABCp1D_Km2 = factory.make_child(a, b, c + _1, d, 0u, oKm2);
          auto ABCD_Km2 = factory.make_child(a, b, c, d, 0u, oKm2);
          if (is_simple()) {
            expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_0") *
                     (Ap1BCD_Km2 - ABCp1D_Km2 +
                      Vector("R12kG12_pfac4_0")[dir] * ABCD_Km2);
            nflops_ += 6;
          }
        }
        return;
      }
      // Build on C
      if (part == 1 && where == InBra) {
        F a(Tint->bra(0, 0));
        F b(Tint->ket(0, 0));
        F c(Tint->bra(1, 0) - _1);
        if (!exists(c)) return;
        F d(Tint->ket(1, 0));

        auto ABCD_K = factory.make_child(a, b, c, d, 0u, oK);
        if (is_simple()) {
          expr_ = Vector("R12kG12_pfac0_1")[dir] * ABCD_K;
          nflops_ += 1;
        }
        const F& cm1 = c - _1;
        if (exists(cm1)) {
          auto ABCm1D_K = factory.make_child(a, b, cm1, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(c[dir]) * Scalar("R12kG12_pfac1_1") * ABCm1D_K;
            nflops_ += 3;
          }
        }
        const F& am1 = a - _1;
        if (exists(am1)) {
          auto Am1BCD_K = factory.make_child(am1, b, c, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(a[dir]) * Scalar("R12kG12_pfac2") * Am1BCD_K;
            nflops_ += 3;
          }
        }
        if (K != 0) {
          const R12kG12 oKm2((R12kG12Descr(K - 2)));
          auto ABCp1D_Km2 = factory.make_child(a, b, c + _1, d, 0u, oKm2);
          auto Ap1BCD_Km2 = factory.make_child(a + _1, b, c, d, 0u, oKm2);
          auto ABCD_Km2 = factory.make_child(a, b, c, d, 0u, oKm2);
          if (is_simple()) {
            expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_1") *
                     (ABCp1D_Km2 - Ap1BCD_Km2 +
                      Vector("R12kG12_pfac4_1")[dir] * ABCD_Km2);
            nflops_ += 6;
          }
        }
        return;
      }
    }  // end of a0c0 case

    if (sxsx) {
      // Build on B
      if (part == 0 && where == InKet) {
        F a(Tint->bra(0, 0));
        F b(Tint->ket(0, 0) - _1);
        if (!exists(b)) return;
        F c(Tint->bra(1, 0));
        F d(Tint->ket(1, 0));

        auto ABCD_K = factory.make_child(a, b, c, d, 0u, oK);
        if (is_simple()) {
          expr_ = Vector("R12kG12_pfac0_0")[dir] * ABCD_K;
          nflops_ += 1;
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          auto ABm1CD_K = factory.make_child(a, bm1, c, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(b[dir]) * Scalar("R12kG12_pfac1_0") * ABm1CD_K;
            nflops_ += 3;
          }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          auto ABCDm1_K = factory.make_child(a, b, c, dm1, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(d[dir]) * Scalar("R12kG12_pfac2") * ABCDm1_K;
            nflops_ += 3;
          }
        }
        if (K != 0) {
          const R12kG12 oKm2((R12kG12Descr(K - 2)));
          auto ABp1CD_Km2 = factory.make_child(a, b + _1, c, d, 0u, oKm2);
          auto ABCDp1_Km2 = factory.make_child(a, b, c, d + _1, 0u, oKm2);
          auto ABCD_Km2 = factory.make_child(a, b, c, d, 0u, oKm2);
          if (is_simple()) {
            expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_0") *
                     (ABp1CD_Km2 - ABCDp1_Km2 +
                      Vector("R12kG12_pfac4_0")[dir] * ABCD_Km2);
            nflops_ += 6;
          }
        }
        return;
      }
      // Build on D
      if (part == 1 && where == InKet) {
        F a(Tint->bra(0, 0));
        F b(Tint->ket(0, 0));
        F c(Tint->bra(1, 0));
        F d(Tint->ket(1, 0) - _1);
        if (!exists(d)) return;

        auto ABCD_K = factory.make_child(a, b, c, d, 0u, oK);
        if (is_simple()) {
          expr_ = Vector("R12kG12_pfac0_1")[dir] * ABCD_K;
          nflops_ += 1;
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          auto ABCDm1_K = factory.make_child(a, b, c, dm1, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(d[dir]) * Scalar("R12kG12_pfac1_1") * ABCDm1_K;
            nflops_ += 3;
          }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          auto ABm1CD_K = factory.make_child(a, bm1, c, d, 0u, oK);
          if (is_simple()) {
            expr_ += Scalar(b[dir]) * Scalar("R12kG12_pfac2") * ABm1CD_K;
            nflops_ += 3;
          }
        }
        if (K != 0) {
          const R12kG12 oKm2((R12kG12Descr(K - 2)));
          auto ABCDp1_Km2 = factory.make_child(a, b, c, d + _1, 0u, oKm2);
          auto ABp1CD_Km2 = factory.make_child(a, b + _1, c, d, 0u, oKm2);
          auto ABCD_Km2 = factory.make_child(a, b, c, d, 0u, oKm2);
          if (is_simple()) {
            expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_1") *
                     (ABCDp1_Km2 - ABp1CD_Km2 +
                      Vector("R12kG12_pfac4_1")[dir] * ABCD_Km2);
            nflops_ += 6;
          }
        }
        return;
      }
    }  // end of 0b0d case

    return;
  }  // K >= 0
}

#if LIBINT_ENABLE_GENERIC_CODE
template <class F, int part, FunctionPosition where>
bool VRR_11_R12kG12_11<F, part, where>::has_generic(
    const std::shared_ptr<CompilationParameters>& cparams) const {
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const unsigned int max_opt_am = cparams->max_am_opt();
  // generic code works for a0c0 and 0a0c classes where am(a) > 1 and am(c) > 1
  // to generate optimized code for xxxx integral need to generate specialized
  // code for up to (x+x)0(x+x)0 integrals
  if (!TrivialBFSet<F>::result && sh_b.zero() && sh_d.zero() &&
      (sh_a.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_c.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_a.norm() > 1u && sh_c.norm() > 1u))
    return true;
  if (!TrivialBFSet<F>::result && sh_a.zero() && sh_c.zero() &&
      (sh_b.norm() > std::max(2 * max_opt_am, 1u) ||
       sh_d.norm() > std::max(2 * max_opt_am, 1u)) &&
      (sh_b.norm() > 1u && sh_d.norm() > 1u))
    return true;
  return false;
}

template <class F, int part, FunctionPosition where>
std::string VRR_11_R12kG12_11<F, part, where>::generic_header() const {
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();
  const int K = target_->oper()->descr().K();
  if (K == -1) {
    if (xsxs) return std::string("OSVRR_xs_xs.h");
    if (sxsx) return std::string("OSVRR_sx_sx.h");
  } else {
    return std::string("VRR_r12kg12_xs_xs.h");
  }
  assert(false);  // unreachable
}

template <class F, int part, FunctionPosition where>
std::string VRR_11_R12kG12_11<F, part, where>::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  const int K = target_->oper()->descr().K();
  std::ostringstream oss;
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const bool xsxs = sh_b.zero() && sh_d.zero();
  const bool sxsx = sh_a.zero() && sh_c.zero();

  bool ahlrichs_simplification = false;
  bool unit_s = false;
  if (xsxs) {
    unit_s = sh_b.is_unit();
  } else {  // (sxsx)
    unit_s = sh_a.is_unit();
  }

  oss << "using namespace libint2;" << std::endl;
  if (K == -1) {
    if (xsxs) {
      oss << "libint2::OSVRR_xs_xs<" << part << "," << sh_a.norm() << ","
          << sh_c.norm() << ",";
      if (not ahlrichs_simplification) oss << (unit_s ? "true," : "false,");
      oss << ((context->cparams()->max_vector_length() == 1) ? "false"
                                                             : "true");
    }
    if (sxsx) {
      oss << "libint2::OSVRR_sx_sx<" << part << "," << sh_b.norm() << ","
          << sh_d.norm() << ",";
      if (not ahlrichs_simplification) oss << (unit_s ? "true," : "false,");
      oss << ((context->cparams()->max_vector_length() == 1) ? "false"
                                                             : "true");
    }
  } else {
    if (xsxs) {
      oss << "libint2::VRR_r12kg12_xs_xs<" << part << "," << sh_a.norm() << ","
          << sh_c.norm() << ",";
      oss << K << ","
          << ((context->cparams()->max_vector_length() == 1) ? "false"
                                                             : "true");
    }
    if (sxsx) {
      // NOTE that using same function, the only difference is in the RR
      // prefactors
      oss << "libint2::VRR_r12kg12_xs_xs<" << part << "," << sh_b.norm() << ","
          << sh_d.norm() << ",";
      oss << K << ","
          << ((context->cparams()->max_vector_length() == 1) ? "false"
                                                             : "true");
    }
  }
  oss << ">::compute(inteval";

  // if K == -1 follow the same logic as for standard VRR
  if (K == -1) {
    oss << "," << args->symbol(0);  // target
    if (not ahlrichs_simplification &&
        unit_s)     // purely to avoid having a 4-term generic RR, reuse 5-term
                    // with dummy argument
      oss << ",0";  // src0-> 0x0
    const unsigned int nargs = args->n();
    for (unsigned int a = 1; a < nargs; a++) {
      oss << "," << args->symbol(a);
    }
  } else {  // K != -1
    const unsigned int nargs = args->n();
    for (unsigned int a = 0; a < nargs; a++) {
      oss << "," << args->symbol(a);
    }
  }
  // if K == 0 add dummy arguments so that the same generic function can be used
  // for all K>=0 cases
  if (K == 0) {
    for (unsigned int a = 0; a < 3; ++a) {
      oss << ",(LIBINT2_REALTYPE*)0";
    }
  }

  oss << ");";

  return oss.str();
}
#endif  // #if !LIBINT_ENABLE_GENERIC_CODE

};  // namespace libint2

#endif
