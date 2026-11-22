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

#include <functional>
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
 *     (\vec{r} - \vec{A})\, \phi_a + (\vec{R}_{AB} \times \vec{B})\, \phi_a
 * \f$, where \f$ (\vec{r} - \vec{A})_j \phi_a = \phi_{a+1_j} \f$, and \f$
 * (\vec{r} - \vec{B})_j \partial_i \phi_b = - \partial \phi_{b+1_j} / \partial
 * B_i \f$ for \f$ i \neq j \f$ .
 *
 * Geometric derivatives (derivative quanta of the Gaussians) are the
 * derivatives of \f$ X_k \f$, including those of the London phase factors,
 * i.e. of \f$ \vec{R}_{AB} \f$ (Leibniz rule), with the identities above
 * generalized to differentiated Gaussians, e.g.
 * \f$ (\vec{r} - \vec{A})_j \partial^n \phi_a / \partial A^n =
 * \partial^n \phi_{a+1_j} / \partial A^n +
 * n_j \partial^{n-1_j} \phi_a / \partial A^{n-1_j} \f$.
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
  /// N.B. this is only a hint (the number of children depends on the
  /// derivative quanta and can exceed it) as long as the result of
  /// ChildFactory::make_child() is copied: it is a reference to an element of
  /// the vector of children, which reallocates beyond this size
  static const unsigned int max_nchildren = 20;

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
  /// @note for shells \p xyz is ignored: a shell contains all components, and
  /// for shell sets (`!is_simple()`) this CR only establishes the dependence on
  /// the child shell sets; the expressions are built when the CR is applied to
  /// the individual integrals (CGF), for which the direction matters
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
  /// @return \p f with one less derivative quantum in direction \p xyz
  static F lower(const F& f, unsigned int xyz) {
    F result(f);
    assert(f.deriv().d(xyz) > 0);
    result.deriv().dec(xyz);
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

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  // (k,l,m) is a cyclic permutation of (x,y,z): (u x v)_k = u_l v_m - u_m v_l
  const auto k = oper->descr().axis();
  const auto l = (k + 1) % 3;
  const auto m = (k + 2) % 3;

  // geometric derivatives: a and b may carry derivative quanta n and p, i.e.
  // the target is d^n/dA^n d^p/dB^p X_k , where the derivatives act on the
  // Gaussians and on R_AB = A - B (the London phase factors); the latter
  // contributes (Leibniz rule; R_AB is linear in A and B)
  //   d^n/dA^n d^p/dB^p (R_AB)_j W = (R_AB)_j W[n,p] + n_j W[n-1_j,p]
  //                                  - p_j W[n,p-1_j]
  // where W[n,p] is W over the differentiated Gaussians; the identities
  // used to remove r from the integrand hold for differentiated Gaussians if
  //   (r-A)_j d^n phi_a/dA^n = d^n phi_{a+1_j}/dA^n + n_j d^{n-1_j} phi_a/dA^n
  const auto& n = a.deriv();
  const auto& p = b.deriv();

  typedef CartesianMultipoleOper<3u> CMultipoleOper;
  ChildFactory<ThisType, GenIntegralSet_1_1<F, CMultipoleOper, EmptySet>>
      cmfactory(this);
  const auto overlap = CMultipoleOper();

  // (a| (r-O)_j |b)
  auto dipole = [&](const F& a, const F& b, unsigned int j) {
    CartesianMultipole_Descr<3u> descr;
    descr.inc(j);
    return cmfactory.make_child(a, b, EmptySet(), CMultipoleOper(descr));
  };
  // (a| dB_i b+1_j) = - (a| (r-B)_j \partial_i |b) for i != j, where dB_i is
  // the derivative w.r.t. B_i (this holds for differentiated b too)
  auto dB_i_bp1_j = [&](unsigned int i, unsigned int j) {
    return cmfactory.make_child(a, differentiate(raise(b, j), i), EmptySet(),
                                overlap);
  };

  if (std::is_same<Oper, GIAOOverlapOper>::value) {
    // X_k = 1/2 [ AB_l (a|(r-O)_m|b) - AB_m (a|(r-O)_l|b) ] + derivatives of AB
    auto r_l = dipole(a, b, l);
    auto r_m = dipole(a, b, m);
    if (is_simple()) {
      expr_ = Scalar(0.5) * (Vector("AB")[l] * r_m - Vector("AB")[m] * r_l);
      nflops_ += 4;
    }
    // + 1/2 [ n_l W_m[n-1_l,p] - n_m W_l[n-1_m,p]
    //         - p_l W_m[n,p-1_l] + p_m W_l[n,p-1_m] ], W_j = (a|(r-O)_j|b)
    struct Term {
      unsigned int count;  // n_i or p_i
      bool on_a;           // differentiated w.r.t. A?
      unsigned int i, j;   // W_j[... - 1_i ...]
      int sign;
    };
    for (const auto& t :
         {Term{n.d(l), true, l, m, +1}, Term{n.d(m), true, m, l, -1},
          Term{p.d(l), false, l, m, -1}, Term{p.d(m), false, m, l, +1}}) {
      if (t.count == 0) continue;
      auto W_j = t.on_a ? dipole(lower(a, t.i), b, t.j)
                        : dipole(a, lower(b, t.i), t.j);
      if (is_simple()) {
        expr_ += Scalar(0.5 * t.sign * t.count) * W_j;
        nflops_ += 2;
      }
    }
  } else if (std::is_same<Oper, GIAOKineticOper>::value ||
             std::is_same<Oper, GIAOElecPotOper>::value) {
    // (a'|O|b') for Gaussians a', b'
    std::function<std::shared_ptr<DGVertex>(const F&, const F&)> O_child;
    ChildFactory<ThisType, GenIntegralSet_1_1<F, KineticOper, EmptySet>>
        tfactory(this);
    ChildFactory<ThisType, GenIntegralSet_1_1<F, ElecPotOper, mType>> vfactory(
        this);
    if (std::is_same<Oper, GIAOKineticOper>::value)
      O_child = [&](const F& a, const F& b) {
        return tfactory.make_child(a, b);
      };
    else
      O_child = [&](const F& a, const F& b) {
        return vfactory.make_child(a, b, mType(0u));
      };

    // Y_j(a',b') = (a'+1_j|O|b') + n'_j (a'-1_j'|O|b') + BO_j (a'|O|b'),
    // where n' are the derivative quanta of a' and a'-1_j' is a' with one
    // less derivative quantum along j, so that
    // (a'|r_j O|b') = Y_j(a',b') + AB_j (a'|O|b') =: W_j(a',b')
    // and AB x W(a,b) = AB x Y(a,b)
    // N.B. returns nullptr if !is_simple(), i.e. only makes the children
    auto Y = [&](const F& a, const F& b, unsigned int j,
                 bool add_AB) -> std::shared_ptr<DGVertex> {
      auto a1j_O_b = O_child(raise(a, j), b);
      auto a_O_b = O_child(a, b);
      std::shared_ptr<DGVertex> am1j_O_b;
      if (a.deriv().d(j) != 0) am1j_O_b = O_child(lower(a, j), b);
      if (!is_simple()) return nullptr;
      std::shared_ptr<RecurrenceRelation::ExprType> result =
          a1j_O_b + Vector("BO")[j] * a_O_b;
      nflops_ += 2;
      if (add_AB) {
        result += Vector("AB")[j] * a_O_b;
        nflops_ += 2;
      }
      if (am1j_O_b) {
        result += Scalar(static_cast<double>(a.deriv().d(j))) * am1j_O_b;
        nflops_ += 2;
      }
      return result;
    };
    // X_k = 1/2 [ AB_l Y_m - AB_m Y_l ] - 1/2 (a|((r-B) x nabla)_k|b)
    auto Y_l = Y(a, b, l, false);
    auto Y_m = Y(a, b, m, false);
    if (is_simple()) {
      expr_ = Scalar(0.5) * (Vector("AB")[l] * Y_m - Vector("AB")[m] * Y_l);
      nflops_ += 4;
    }
    // + 1/2 [ n_l W_m[n-1_l,p] - n_m W_l[n-1_m,p]
    //         - p_l W_m[n,p-1_l] + p_m W_l[n,p-1_m] ]
    struct Term {
      unsigned int count;  // n_i or p_i
      bool on_a;           // differentiated w.r.t. A?
      unsigned int i, j;   // W_j[... - 1_i ...]
      int sign;
    };
    for (const auto& t :
         {Term{n.d(l), true, l, m, +1}, Term{n.d(m), true, m, l, -1},
          Term{p.d(l), false, l, m, -1}, Term{p.d(m), false, m, l, +1}}) {
      if (t.count == 0) continue;
      auto W_j = t.on_a ? Y(lower(a, t.i), b, t.j, true)
                        : Y(a, lower(b, t.i), t.j, true);
      if (is_simple()) {
        expr_ += Scalar(0.5 * t.sign * t.count) * W_j;
        nflops_ += 2;
      }
    }

    if (std::is_same<Oper, GIAOKineticOper>::value) {
      // (a|((r-B) x nabla)_k|b) = (a|(r-B)_l \partial_m|b)
      //                           - (a|(r-B)_m \partial_l|b)
      //                         = (a|dB_l b+1_m) - (a|dB_m b+1_l)
      // the dependence of (r-B) on B is accounted for by differentiating
      // b+1_j (the identity (r-B)_j phi_b = phi_{b+1_j} holds for any B)
      auto dl_b1m = dB_i_bp1_j(l, m);
      auto dm_b1l = dB_i_bp1_j(m, l);
      if (is_simple()) {
        expr_ -= Scalar(0.5) * (dl_b1m - dm_b1l);
        nflops_ += 3;
      }
    }
  } else if (std::is_same<Oper, AngularMomentumOper>::value) {
    // (r-O)_j d_i b = (r-B)_j d_i b + BO_j d_i b, where d_i b = - dB_i b, hence
    // X_k = (a|(r-O)_l d_m b) - (a|(r-O)_m d_l b)
    //     = (a|dB_l b+1_m) - (a|dB_m b+1_l)
    //       + BO_m (a|dB_l b) - BO_l (a|dB_m b)
    //       + p_m (a|dB_l (b-1_m')) - p_l (a|dB_m (b-1_l'))
    // the last line is due to the derivative quanta p of b (the integrand
    // depends on A and B only via the Gaussians, O is fixed), see
    // (r-B)_j d^p phi_b/dB^p above
    auto dl_b1m = dB_i_bp1_j(l, m);
    auto dm_b1l = dB_i_bp1_j(m, l);
    auto dl_b =
        cmfactory.make_child(a, differentiate(b, l), EmptySet(), overlap);
    auto dm_b =
        cmfactory.make_child(a, differentiate(b, m), EmptySet(), overlap);
    if (is_simple()) {
      expr_ =
          (dl_b1m - dm_b1l) + (Vector("BO")[m] * dl_b - Vector("BO")[l] * dm_b);
      nflops_ += 5;
    }
    if (p.d(m) != 0) {
      auto dl_bm1m = cmfactory.make_child(a, differentiate(lower(b, m), l),
                                          EmptySet(), overlap);
      if (is_simple()) {
        expr_ += Scalar(static_cast<double>(p.d(m))) * dl_bm1m;
        nflops_ += 2;
      }
    }
    if (p.d(l) != 0) {
      auto dm_bm1l = cmfactory.make_child(a, differentiate(lower(b, l), m),
                                          EmptySet(), overlap);
      if (is_simple()) {
        expr_ -= Scalar(static_cast<double>(p.d(l))) * dm_bm1l;
        nflops_ += 2;
      }
    }
  } else
    throw std::logic_error("CR_1_GIAO_1: unsupported operator");

}  // CR_1_GIAO_1<F,Oper>::CR_1_GIAO_1

};  // namespace libint2

#endif  // LIBINT_COMP_1_GIAO_1_H
