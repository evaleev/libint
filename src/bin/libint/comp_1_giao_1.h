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

#ifndef LIBINT_COMP_1_GIAO_1_H
#define LIBINT_COMP_1_GIAO_1_H

#include <generic_rr.h>

#include <type_traits>

namespace libint2 {

/**
 * Computes the magnetic one-body integrals over Cartesian Gaussians by
 * rewriting them as linear combinations of conventional integrals (overlap,
 * multipole, kinetic, electrostatic potential) over Gaussians with raised
 * angular momentum and/or differentiated with respect to their centers.
 *
 * Zero-field derivatives of integrals over London orbitals (GIAOs)
 * \f$ \omega_a = e^{-\frac{i}{2} (\vec{B} \times \vec{A}) \cdot \vec{r}}
 * \phi_a \f$ with respect to the magnetic field \f$ \vec{B} \f$ are purely
 * imaginary: \f$ \partial O_{ab} / \partial B_k = i X_k \f$ (Helgaker and
 * Jørgensen, J. Chem. Phys. 95, 2595 (1991)). The real \f$ X_k \f$ are computed
 * using
 * \f$ (\vec{R}_{AB} \times \vec{r})\, \phi_a = \vec{R}_{AB} \times
 *     (\vec{r} - \vec{A})\, \phi_a + (\vec{R}_{AB} \times \vec{B})\, \phi_a \f$,
 * where \f$ (\vec{r} - \vec{A})_j \phi_a = \phi_{a+1_j} \f$, and
 * \f$ (\vec{r} - \vec{B})_j \partial_i \phi_b = - \partial \phi_{b+1_j} /
 * \partial B_i \f$ for \f$ i \neq j \f$ .
 *
 * The \c BO prefactor (\f$ \vec{B} - \vec{O} \f$) is the position of the ket
 * center relative to the origin (the coordinate origin for the GIAO
 * operators, the gauge origin for AngularMomentumOper).
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 * @tparam Oper one of GIAOOverlapOper, GIAOKineticOper, GIAOElecPotOper, or
 * AngularMomentumOper
 */
template <typename F, typename Oper>
class CR_1_GIAO_1
    : public GenericRecurrenceRelation<CR_1_GIAO_1<F, Oper>, F,
                                       GenIntegralSet_1_1<F, Oper, EmptySet>> {
 public:
  typedef CR_1_GIAO_1<F, Oper> ThisType;
  typedef F BasisFunctionType;
  typedef Oper OperType;
  typedef GenIntegralSet_1_1<F, Oper, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 8;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_1_GIAO_1(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  /// @return \p f with its angular momentum raised by 1 in direction \p xyz
  /// (for shells: angular momentum raised by 1)
  static F raise(const F& f, unsigned int xyz) {
    F result(f);
    result.inc(std::is_same<F, CGShell>::value ? 0 : xyz);
    return result;
  }
  /// @return \p f differentiated with respect to its center in direction
  /// \p xyz
  static F differentiate(const F& f, unsigned int xyz) {
    F result(f);
    result.deriv().inc(xyz);
    return result;
  }
};

template <typename F, typename Oper>
CR_1_GIAO_1<F, Oper>::CR_1_GIAO_1(const std::shared_ptr<TargetType>& Tint,
                                  unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  const auto& a = Tint->bra(0, 0);
  const auto& b = Tint->ket(0, 0);
  const auto& oper = Tint->oper();

  // like CR_1_σpVσp_1, only primitive integrals are rewritten; contracted
  // integrals are uncontracted first
  if (a.contracted() || b.contracted()) return;

  // prefactors that depend on the Gaussian centers are not differentiated, so
  // geometric derivatives of these integrals cannot be expressed this way
  if (a.deriv().norm() != 0 || b.deriv().norm() != 0)
    throw std::logic_error(
        "CR_1_GIAO_1: geometric derivatives of GIAO/angular momentum integrals "
        "are not supported");

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  // (k,l,m) is a cyclic permutation of (x,y,z): (u x v)_k = u_l v_m - u_m v_l
  const auto k = oper->descr().axis();
  const auto l = (k + 1) % 3;
  const auto m = (k + 2) % 3;

  typedef CartesianMultipoleOper<3u> CMultipoleOper;
  ChildFactory<ThisType, GenIntegralSet_1_1<F, CMultipoleOper, EmptySet>>
      cmfactory(this);
  const auto overlap = CMultipoleOper();

  // (a| (r-O)_j |b)
  auto dipole = [&](unsigned int j) {
    CartesianMultipole_Descr<3u> descr;
    descr.inc(j);
    return cmfactory.make_child(a, b, EmptySet(), CMultipoleOper(descr));
  };
  // (a| dB_i b+1_j) = - (a| (r-B)_j \partial_i |b) for i != j, where dB_i is
  // the derivative w.r.t. B_i
  auto dB_i_bp1_j = [&](unsigned int i, unsigned int j) {
    return cmfactory.make_child(a, differentiate(raise(b, j), i), EmptySet(),
                                overlap);
  };

  if (std::is_same<Oper, GIAOOverlapOper>::value) {
    // X_k = 1/2 [ AB_l (a|(r-O)_m|b) - AB_m (a|(r-O)_l|b) ]
    auto r_l = dipole(l);
    auto r_m = dipole(m);
    if (is_simple()) {
      expr_ = Scalar(0.5) *
              (Vector("AB")[l] * r_m - Vector("AB")[m] * r_l);
      nflops_ += 4;
    }
  } else if (std::is_same<Oper, GIAOKineticOper>::value ||
             std::is_same<Oper, GIAOElecPotOper>::value) {
    // Y_j = (a+1_j|O|b) + BO_j (a|O|b), so that (a|((r-O) O)_j|b) = Y_j +
    // AB_j (a|O|b) and AB x (a|(r-O) O|b) = AB x Y
    std::shared_ptr<DGVertex> a1l_O_b, a1m_O_b, a_O_b;
    if (std::is_same<Oper, GIAOKineticOper>::value) {
      ChildFactory<ThisType, GenIntegralSet_1_1<F, KineticOper, EmptySet>>
          factory(this);
      a1l_O_b = factory.make_child(raise(a, l), b);
      a1m_O_b = factory.make_child(raise(a, m), b);
      a_O_b = factory.make_child(a, b);
    } else {
      ChildFactory<ThisType, GenIntegralSet_1_1<F, ElecPotOper, mType>>
          factory(this);
      const mType zero_m(0u);
      a1l_O_b = factory.make_child(raise(a, l), b, zero_m);
      a1m_O_b = factory.make_child(raise(a, m), b, zero_m);
      a_O_b = factory.make_child(a, b, zero_m);
    }

    // X_k = 1/2 [ AB_l Y_m - AB_m Y_l ] - 1/2 (a|((r-B) x nabla)_k|b)
    // (a|((r-B) x nabla)_k|b) = (a|(r-B)_l \partial_m|b)
    //                           - (a|(r-B)_m \partial_l|b)
    //                         = (a|dB_l b+1_m) - (a|dB_m b+1_l)
    std::shared_ptr<DGVertex> dl_b1m, dm_b1l;
    if (std::is_same<Oper, GIAOKineticOper>::value) {
      dl_b1m = dB_i_bp1_j(l, m);
      dm_b1l = dB_i_bp1_j(m, l);
    }
    if (is_simple()) {
      auto Y_l = a1l_O_b + Vector("BO")[l] * a_O_b;
      auto Y_m = a1m_O_b + Vector("BO")[m] * a_O_b;
      auto AB_x_Y = Vector("AB")[l] * Y_m - Vector("AB")[m] * Y_l;
      nflops_ += 7;
      if (std::is_same<Oper, GIAOKineticOper>::value) {
        expr_ = Scalar(0.5) * (AB_x_Y - (dl_b1m - dm_b1l));
        nflops_ += 3;
      } else {
        expr_ = Scalar(0.5) * AB_x_Y;
        nflops_ += 1;
      }
    }
  } else if (std::is_same<Oper, AngularMomentumOper>::value) {
    // (r-O)_j d_i b = (r-B)_j d_i b + BO_j d_i b, where d_i b = - dB_i b, hence
    // X_k = (a|(r-O)_l d_m b) - (a|(r-O)_m d_l b)
    //     = (a|dB_l b+1_m) - (a|dB_m b+1_l)
    //       + BO_m (a|dB_l b) - BO_l (a|dB_m b)
    auto dl_b1m = dB_i_bp1_j(l, m);
    auto dm_b1l = dB_i_bp1_j(m, l);
    auto dl_b = cmfactory.make_child(a, differentiate(b, l), EmptySet(),
                                     overlap);
    auto dm_b = cmfactory.make_child(a, differentiate(b, m), EmptySet(),
                                     overlap);
    if (is_simple()) {
      expr_ = (dl_b1m - dm_b1l) +
              (Vector("BO")[m] * dl_b - Vector("BO")[l] * dm_b);
      nflops_ += 5;
    }
  } else
    throw std::logic_error("CR_1_GIAO_1: unsupported operator");

}  // CR_1_GIAO_1<F,Oper>::CR_1_GIAO_1

};  // namespace libint2

#endif  // LIBINT_COMP_1_GIAO_1_H
