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

#ifndef LIBINT_COMP_11_GIAO_11_H
#define LIBINT_COMP_11_GIAO_11_H

#include <generic_rr.h>

#include <type_traits>

namespace libint2 {

/**
 * Computes the zero-field derivatives of the electron repulsion integrals over
 * London orbitals (GIAOs) with respect to the magnetic field (see
 * GIAOTwoPRep_Descr) as linear combinations of conventional electron repulsion
 * integrals with raised angular momentum.
 *
 * The London phase factors of the bra pair (ab| of electron 1 combine into
 * \f$ e^{\frac{i}{2} \vec{B} \cdot (\vec{R}_{AB} \times \vec{r}_1)} \f$ (and
 * similarly for the ket pair |cd) of electron 2), hence
 * \f$ X^{(1)}_k = \frac{1}{2} (\vec{R}_{AB} \times \vec{Y})_k \f$,
 * \f$ Y_j = (a+1_j\, b|cd) + B_j (ab|cd) \f$, where
 * \f$ \vec{R}_{AB} \times \vec{r} = \vec{R}_{AB} \times (\vec{r} - \vec{A}) +
 * \vec{R}_{AB} \times \vec{B} \f$ was used; \f$ X^{(2)}_k \f$ is obtained
 * analogously with \f$ (ab|c+1_j\, d) \f$ and \f$ \vec{D} \f$. The positions
 * of the second center of each pair, \f$ \vec{B} \f$ and \f$ \vec{D} \f$,
 * relative to the coordinate origin are the \c BO and \c DO prefactors.
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_GIAOTwoPRep_11
    : public GenericRecurrenceRelation<
          CR_11_GIAOTwoPRep_11<F>, F,
          GenIntegralSet_11_11<F, GIAOTwoPRep, mType>> {
 public:
  typedef CR_11_GIAOTwoPRep_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef GIAOTwoPRep OperType;
  typedef GenIntegralSet_11_11<F, GIAOTwoPRep, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 3;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_11_GIAOTwoPRep_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  /// @return \p f with its angular momentum raised by 1 in direction \p xyz
  /// (for shells: angular momentum raised by 1)
  static F raise(const F& f, unsigned int xyz) {
    F result(f);
    result.inc(std::is_same<F, CGShell>::value ? 0 : xyz);
    return result;
  }
};

template <typename F>
CR_11_GIAOTwoPRep_11<F>::CR_11_GIAOTwoPRep_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int)
    : ParentType(Tint, 0) {
  const auto& a = Tint->bra(0, 0);
  const auto& b = Tint->ket(0, 0);
  const auto& c = Tint->bra(1, 0);
  const auto& d = Tint->ket(1, 0);
  const auto& aux = Tint->aux();
  const auto& oper = Tint->oper();

  // prefactors that depend on the Gaussian centers are not differentiated, so
  // geometric derivatives of these integrals cannot be expressed this way
  if (a.deriv().norm() != 0 || b.deriv().norm() != 0 ||
      c.deriv().norm() != 0 || d.deriv().norm() != 0)
    throw std::logic_error(
        "CR_11_GIAOTwoPRep_11: geometric derivatives of GIAO electron "
        "repulsion integrals are not supported");

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  // (k,l,m) is a cyclic permutation of (x,y,z): (u x v)_k = u_l v_m - u_m v_l
  const auto k = oper->descr().axis();
  const auto l = (k + 1) % 3;
  const auto m = (k + 2) % 3;
  const auto bra_pair = oper->descr().particle() == 0;

  ChildFactory<ThisType, GenIntegralSet_11_11<F, TwoPRep, mType>> factory(
      this);
  auto abcd = factory.make_child(a, b, c, d, aux);
  // (a+1_j b|cd) or (ab|c+1_j d)
  auto raised = [&](unsigned int j) {
    return bra_pair ? factory.make_child(raise(a, j), b, c, d, aux)
                    : factory.make_child(a, b, raise(c, j), d, aux);
  };
  auto raised_l = raised(l);
  auto raised_m = raised(m);

  if (is_simple()) {
    // R = AB, O = BO for the bra pair, R = CD, O = DO for the ket pair
    auto R = bra_pair ? Vector("AB") : Vector("CD");
    auto O = bra_pair ? Vector("BO") : Vector("DO");
    auto Y_l = raised_l + O[l] * abcd;
    auto Y_m = raised_m + O[m] * abcd;
    expr_ = Scalar(0.5) * (R[l] * Y_m - R[m] * Y_l);
    nflops_ += 8;
  }
}  // CR_11_GIAOTwoPRep_11<F>::CR_11_GIAOTwoPRep_11

};  // namespace libint2

#endif  // LIBINT_COMP_11_GIAO_11_H
