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

#ifndef _libint2_src_bin_libint_compxyz_h_
#define _libint2_src_bin_libint_compxyz_h_

#include <generic_rr.h>

namespace libint2 {

/**
 * this computes integral of @c Oper over CGShell/CGF as a product of 1-d
 * integrals
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F, typename Oper, typename AuxQuanta = EmptySet>
class CR_XYZ_1_1 : public GenericRecurrenceRelation<
                       CR_XYZ_1_1<F, Oper, AuxQuanta>, F,
                       GenIntegralSet_1_1<F, Oper, AuxQuanta> > {
 public:
  typedef CR_XYZ_1_1<F, Oper, AuxQuanta> ThisType;
  typedef F BasisFunctionType;
  typedef Oper OperType;
  typedef GenIntegralSet_1_1<F, Oper, AuxQuanta> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 100;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_XYZ_1_1(const std::shared_ptr<TargetType>&, unsigned int dir = 0);

  static std::string descr() { return "CR"; }

  /// specialize this for the given operator type and CGF
  void compute(const BasisFunctionType& bra, const BasisFunctionType& ket,
               const Oper& oper);
};

template <typename F, typename Oper, typename AuxQuanta>
CR_XYZ_1_1<F, Oper, AuxQuanta>::CR_XYZ_1_1(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  // WARNING assuming one function per position
  const auto& a = Tint->bra(0, 0);
  const auto& b = Tint->ket(0, 0);
  const auto& aux = Tint->aux();
  const auto& oper = Tint->oper();

  {
    // can't apply to contracted basis functions
    if (a.contracted() || b.contracted()) return;
    // can't apply to differentiated CGF (derivatives will be expanded first)
    if (TrivialBFSet<F>::result &&
        (a.deriv().norm() != 0 || b.deriv().norm() != 0))
      return;
  }

  compute(a, b, oper);
}  // CR_XYZ_1_1<F,Oper,AuxQuanta>::CR_XYZ_1_1

};  // namespace libint2

#endif
