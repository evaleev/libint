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
 * Geometric derivatives (derivative quanta of the Gaussians) are the
 * derivatives of \f$ X_k \f$, including those of the London phase factors,
 * i.e. of \f$ \vec{R}_{AB} \f$ and \f$ \vec{R}_{CD} \f$; they are
 * evaluated as in CR_1_GIAO_1: by the Leibniz rule, with
 * \f$ (\vec{r} - \vec{A})_j \partial^n \phi_a / \partial A^n =
 * \partial^n \phi_{a+1_j} / \partial A^n +
 * n_j \partial^{n-1_j} \phi_a / \partial A^{n-1_j} \f$. Since the phase
 * factors refer to the (fixed) coordinate origin these integrals are not
 * translationally invariant, hence the derivatives with respect to all 4
 * centers are evaluated explicitly (the conventional integrals they are
 * expressed in terms of do use the translational invariance).
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_GIAOTwoPRep_11 : public GenericRecurrenceRelation<
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
  CR_11_GIAOTwoPRep_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

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
  /// @return \p f with one less derivative quantum in direction \p xyz
  static F lower(const F& f, unsigned int xyz) {
    F result(f);
    assert(f.deriv().d(xyz) > 0);
    result.deriv().dec(xyz);
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

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  // (k,l,m) is a cyclic permutation of (x,y,z): (u x v)_k = u_l v_m - u_m v_l
  const auto k = oper->descr().axis();
  const auto l = (k + 1) % 3;
  const auto m = (k + 2) % 3;
  const auto bra_pair = oper->descr().particle() == 0;

  // the pair whose London phase factors this component accounts for consists
  // of functions f and g, on centers F and G: (f,g) = (a,b) for the bra pair,
  // (c,d) for the ket pair; the functions of the other pair are spectators
  // (their derivative quanta, if any, are simply passed on to the children)
  const auto& f = bra_pair ? a : c;
  const auto& g = bra_pair ? b : d;
  // R = F - G, O = G - origin
  auto R = bra_pair ? Vector("AB") : Vector("CD");
  auto O = bra_pair ? Vector("BO") : Vector("DO");

  ChildFactory<ThisType, GenIntegralSet_11_11<F, TwoPRep, mType>> factory(this);
  // the conventional integral with f and g replaced by f' and g'
  auto eri = [&](const F& f, const F& g) {
    return bra_pair ? factory.make_child(f, g, c, d, aux)
                    : factory.make_child(a, b, f, g, aux);
  };

  // geometric derivatives: f and g may carry derivative quanta n and p, i.e.
  // the target is d^n/dF^n d^p/dG^p X_k , where the derivatives act on the
  // Gaussians and on R = F - G (the London phase factors); the latter
  // contributes (Leibniz rule; R is linear in F and G)
  //   d^n/dF^n d^p/dG^p R_j W = R_j W[n,p] + n_j W[n-1_j,p] - p_j W[n,p-1_j]
  // where W[n,p] is W over the differentiated Gaussians, see CR_1_GIAO_1
  const auto& n = f.deriv();
  const auto& p = g.deriv();

  // Y_j(f',g') = (f'+1_j g'|..) + n'_j (f'-1_j' g'|..) + O_j (f' g'|..),
  // where n' are the derivative quanta of f' and f'-1_j' is f' with one less
  // derivative quantum along j, so that
  // (f' r_j g'|..) = Y_j(f',g') + R_j (f' g'|..) =: W_j(f',g')
  // and R x W(f,g) = R x Y(f,g)
  // N.B. returns nullptr if !is_simple(), i.e. only makes the children
  auto Y = [&](const F& f, const F& g, unsigned int j,
               bool add_R) -> std::shared_ptr<RecurrenceRelation::ExprType> {
    auto f1j_g = eri(raise(f, j), g);
    auto f_g = eri(f, g);
    std::shared_ptr<DGVertex> fm1j_g;
    if (f.deriv().d(j) != 0) fm1j_g = eri(lower(f, j), g);
    if (!is_simple()) return nullptr;
    std::shared_ptr<RecurrenceRelation::ExprType> result = f1j_g + O[j] * f_g;
    nflops_ += 2;
    if (add_R) {
      result += R[j] * f_g;
      nflops_ += 2;
    }
    if (fm1j_g) {
      result += Scalar(static_cast<double>(f.deriv().d(j))) * fm1j_g;
      nflops_ += 2;
    }
    return result;
  };

  // X_k = 1/2 [ R_l Y_m - R_m Y_l ]
  auto Y_l = Y(f, g, l, false);
  auto Y_m = Y(f, g, m, false);
  if (is_simple()) {
    expr_ = Scalar(0.5) * (R[l] * Y_m - R[m] * Y_l);
    nflops_ += 4;
  }
  // + 1/2 [ n_l W_m[n-1_l,p] - n_m W_l[n-1_m,p]
  //         - p_l W_m[n,p-1_l] + p_m W_l[n,p-1_m] ]
  struct Term {
    unsigned int count;  // n_i or p_i
    bool on_f;           // differentiated w.r.t. F?
    unsigned int i, j;   // W_j[... - 1_i ...]
    int sign;
  };
  for (const auto& t :
       {Term{n.d(l), true, l, m, +1}, Term{n.d(m), true, m, l, -1},
        Term{p.d(l), false, l, m, -1}, Term{p.d(m), false, m, l, +1}}) {
    if (t.count == 0) continue;
    auto W_j = t.on_f ? Y(lower(f, t.i), g, t.j, true)
                      : Y(f, lower(g, t.i), t.j, true);
    if (is_simple()) {
      expr_ += Scalar(0.5 * t.sign * t.count) * W_j;
      nflops_ += 2;
    }
  }
}  // CR_11_GIAOTwoPRep_11<F>::CR_11_GIAOTwoPRep_11

};  // namespace libint2

#endif  // LIBINT_COMP_11_GIAO_11_H
