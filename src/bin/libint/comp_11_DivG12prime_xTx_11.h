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

#ifndef _libint2_src_bin_libint_cr11divg12primextx11_h_
#define _libint2_src_bin_libint_cr11divg12primextx11_h_

#include <gaussoper.h>
#include <generic_rr.h>
#include <integral_11_11.h>

#define USE_R12kR12lG12 1

namespace libint2 {

/** Compute relation for 2-e integrals of the DivG12prime_xTx operators.
 */
template <class BFSet>
class CR_11_DivG12prime_xTx_11
    : public GenericRecurrenceRelation<
          CR_11_DivG12prime_xTx_11<BFSet>, BFSet,
          GenIntegralSet_11_11<BFSet, DivG12prime_xTx, mType> > {
 public:
  typedef CR_11_DivG12prime_xTx_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, DivG12prime_xTx, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 36;

  using ParentType::Instance;

  /// This relation is not directional
  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_11_DivG12prime_xTx_11(const std::shared_ptr<TargetType>&,
                           unsigned int dir);
  static std::string descr() { return "CR"; }

  template <class RR, class C>
  friend class ChildFactory;
};

template <class F>
CR_11_DivG12prime_xTx_11<F>::CR_11_DivG12prime_xTx_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // kinetic energy of which electron?
  const int i = target_->oper()->descr().I();
  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  if (a.contracted() || b.contracted() || c.contracted() || d.contracted() ||
      target_->oper()->descr().contracted())
    return;

#if USE_R12kR12lG12
  if (i == 0) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
        ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    for (int bxyz = 0; bxyz < 3; ++bxyz) {
      for (int kxyz = 0; kxyz < 3; ++kxyz) {
        R12k_R12l_G12_Descr descr(unit_intvec3(bxyz), unit_intvec3(kxyz));
        factory.wedge(Nabla1(_pbra(a, c), bxyz), Nabla1(_pket(b, d), kxyz),
                      EmptySet(), R12kR12lG12(descr));
      }
    }
  }
  if (i == 1) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kR12lG12, EmptySet>
        ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    for (int bxyz = 0; bxyz < 3; ++bxyz) {
      for (int kxyz = 0; kxyz < 3; ++kxyz) {
        R12k_R12l_G12_Descr descr(unit_intvec3(bxyz), unit_intvec3(kxyz));
        factory.wedge(Nabla2(_pbra(a, c), bxyz), Nabla2(_pket(b, d), kxyz),
                      EmptySet(), R12kR12lG12(descr));
      }
    }
  }
#else
  // |Nabla1.G12' ac ) ^ |Nabla1.G12' bd )
  if (i == 0) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType> ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    factory.wedge(R12vec_dot_Nabla1(_pbra(a, c)),
                  R12vec_dot_Nabla1(_pket(b, d)), mType(0u), R12kG12(0));
  }
  // |Nabla2.G12' ac ) ^ |Nabla2.G12' bd )
  if (i == 1) {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType> ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    factory.wedge(R12vec_dot_Nabla2(_pbra(a, c)),
                  R12vec_dot_Nabla2(_pket(b, d)), mType(0u), R12kG12(0));
  }
#endif
  if (is_simple())
    expr_ *= Scalar(-4.0) * Scalar("gamma_bra") * Scalar("gamma_ket");
}

};  // namespace libint2

#endif
