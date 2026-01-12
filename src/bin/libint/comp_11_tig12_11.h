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

#ifndef _libint2_src_bin_libint_cr11tig1211_h_
#define _libint2_src_bin_libint_cr11tig1211_h_

#include <gaussoper.h>
#include <generic_rr.h>
#include <tig12_11_11.h>

namespace libint2 {

/** Compute relation for 2-e integrals of the Ti_G12 operators.
 */
template <class BFSet>
class CR_11_TiG12_11 : public GenericRecurrenceRelation<
                           CR_11_TiG12_11<BFSet>, BFSet,
                           GenIntegralSet_11_11<BFSet, TiG12, mType> > {
 public:
  typedef CR_11_TiG12_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, TiG12, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 8;

  using ParentType::Instance;

  /// This relation is not directional
  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;
  template <class RR, class C>
  friend class ChildFactory;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_11_TiG12_11(const std::shared_ptr<TargetType>&, unsigned int dir);
  static std::string descr() { return "CR"; }
};

template <class F>
CR_11_TiG12_11<F>::CR_11_TiG12_11(const std::shared_ptr<TargetType>& Tint,
                                  unsigned int dir)
    : ParentType(Tint, dir) {
  if (dir != 0) return;
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  // kinetic energy of which electron?
  const int i = target_->oper()->descr().K();
  const R12kG12 G0(R12_k_G12_Descr(0));
  const R12kG12 G2(R12_k_G12_Descr(2));

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  if (i == 0) {
    if (b.contracted() || target_->oper()->descr().contracted()) return;
  }
  if (i == 1) {
    if (d.contracted() || target_->oper()->descr().contracted()) return;
  }

  // [T1,G12]
  if (i == 0) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
        ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    for (int xyz = 0; xyz < 3; ++xyz) {
      R12k_R12l_G12_Descr descr(IntVec3(), unit_intvec3(xyz));
      factory.wedge(_pbra(a, c), Nabla1(_pket(b, d), xyz), EmptySet(),
                    R12kR12lG12(descr));
    }
    if (is_simple()) expr_ *= Scalar(2.0) * Scalar("gamma");
  }
  // [T2,G12]
  if (i == 1) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
        ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    for (int xyz = 0; xyz < 3; ++xyz) {
      R12k_R12l_G12_Descr descr(IntVec3(), unit_intvec3(xyz));
      factory.wedge(_pbra(a, c), Nabla2(_pket(b, d), xyz), EmptySet(),
                    R12kR12lG12(descr));
    }
    if (is_simple()) expr_ *= Scalar(-2.0) * Scalar("gamma");
  }

  {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType> ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    auto ab_G0_cd = factory.make_child(a, b, c, d, 0u, G0);
    if (is_simple()) expr_ += Scalar(3.0) * Scalar("gamma") * ab_G0_cd;

    auto ab_G2_cd = factory.make_child(a, b, c, d, 0u, G2);
    if (is_simple())
      expr_ += Scalar(-2.0) * Scalar("gamma") * Scalar("gamma") * ab_G2_cd;
  }
}

};  // namespace libint2

#endif
