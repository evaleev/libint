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

#ifndef _libint2_src_bin_libint_cr11r12kr12lg1211_h_
#define _libint2_src_bin_libint_cr11r12kr12lg1211_h_

#include <gaussoper.h>
#include <generic_rr.h>
#include <integral_11_11.h>

namespace libint2 {

/** Compute relation for integrals of operator R12k_R12l_G12.
    Reduce power of R12 applied to bra first, then ket.
    Children are either integrals of R12k_R12l_G12 or R12kG12 with K=0.
*/
template <class BFSet>
class CR_11_R12kR12lG12_11
    : public GenericRecurrenceRelation<
          CR_11_R12kR12lG12_11<BFSet>, BFSet,
          GenIntegralSet_11_11<BFSet, R12kR12lG12, EmptySet> > {
 public:
  typedef CR_11_R12kR12lG12_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, R12kR12lG12, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 3;

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
  CR_11_R12kR12lG12_11(const std::shared_ptr<TargetType>&, unsigned int dir);
  static std::string descr() { return "CR"; }
};

template <class F>
CR_11_R12kR12lG12_11<F>::CR_11_R12kR12lG12_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // std::cout << "CR_11_R12kR12lG12_11<F> -- applying to "
  //           << Tint->label() << std::endl;

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  const IntVec3& r12pbra = Tint->oper()->descr().K();
  const IntVec3& r12pket = Tint->oper()->descr().L();
  const int norm = r12pbra.norm1() + r12pket.norm1();

  // Try reducing r12pbra
  for (int xyz = 0; xyz < 3; ++xyz) {
    const IntVec3& _1 = unit_intvec3(xyz);
    const IntVec3 r12pbra_m1 = r12pbra - _1;

    if (!ltzero(r12pbra_m1)) {
      if (norm == 1) {
        typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType>
            ChildType;
        ChildFactory<ThisType, ChildType> factory(this);
        factory.wedge(R12v(_pbra(a, c), xyz), _pket(b, d), mType(0u),
                      R12kG12(R12_k_G12_Descr(0)));
      } else {
        typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
            ChildType;
        ChildFactory<ThisType, ChildType> factory(this);
        R12k_R12l_G12_Descr odescr(r12pbra_m1, r12pket);
        factory.wedge(R12v(_pbra(a, c), xyz), _pket(b, d), EmptySet(),
                      R12kR12lG12(odescr));
      }
      return;
    }
  }

  // Try reducing r12pket
  for (int xyz = 0; xyz < 3; ++xyz) {
    const IntVec3& _1 = unit_intvec3(xyz);
    const IntVec3 r12pket_m1 = r12pket - _1;

    if (!ltzero(r12pket_m1)) {
      if (norm == 1) {
        typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType>
            ChildType;
        ChildFactory<ThisType, ChildType> factory(this);
        factory.wedge(_pbra(a, c), R12v(_pket(b, d), xyz), mType(0u),
                      R12kG12(R12_k_G12_Descr(0)));
      } else {
        typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
            ChildType;
        ChildFactory<ThisType, ChildType> factory(this);
        R12k_R12l_G12_Descr odescr(r12pbra, r12pket_m1);
        factory.wedge(_pbra(a, c), R12v(_pket(b, d), xyz), EmptySet(),
                      R12kR12lG12(odescr));
      }
      return;
    }
  }
}

};  // namespace libint2

#endif
