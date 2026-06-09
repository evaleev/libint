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

#ifndef LIBINT_RKB_FOLD_1BODY_H
#define LIBINT_RKB_FOLD_1BODY_H

#include <generic_rr.h>

#include <stdexcept>

namespace libint2 {

/**
 * Shared one-body \f$ \sigma\cdot\hat{p}\, \hat{O}\, \sigma\cdot\hat{p} \f$
 * Pauli-quaternion fold for primitive Gaussians. Both \c CR_1_σpVσp_1 (inner
 * operator = electrostatic potential V) and \c CR_1_σpeμσp_1 (inner operator =
 * Cartesian dipole \f$ r_k \f$) reduce to exactly this combination of the 6
 * first-derivative basis functions; only the inner child integral differs, so
 * the caller supplies it via @p make_child.
 *
 * Using \f$ \sigma_a \sigma_b = \delta_{ab} + i\epsilon_{abc}\sigma_c \f$, the
 * quaternion index @p q selects:
 *   - q=0 (trace δ_ab):   Σ_a (∂_a a | O | ∂_a b)
 *   - q=1 (σ_x antisym):  (∂_y a | O | ∂_z b) − (∂_z a | O | ∂_y b)
 *   - q=2 (σ_y antisym):  (∂_z a | O | ∂_x b) − (∂_x a | O | ∂_z b)
 *   - q=3 (σ_z antisym):  (∂_x a | O | ∂_y b) − (∂_y a | O | ∂_x b)
 *
 * Children are always created (so they register on the graph); @p expr is set
 * only when @p is_simple.
 *
 * @tparam F          basis function type (CGShell or CGF)
 * @tparam ExprPtr    type of the RR's expr_ member (shared_ptr<ExprType>)
 * @tparam MakeChild  callable (const F& Da, const F& Db) -> child vertex,
 *                    wrapping the operator-specific ChildFactory::make_child
 * @param a,b         bra/ket basis functions of the target
 * @param q           Pauli quaternion index (0..3)
 * @param is_simple   value of the RR's is_simple()
 * @param expr        the RR's expr_ (assigned when is_simple)
 * @param nflops      the RR's nflops_ (incremented by the fold's adds)
 * @param make_child  builds the inner child integral (Da | O | Db)
 */
template <typename F, typename ExprPtr, typename MakeChild>
void sigma_p_1body_fold(const F &a, const F &b, unsigned int q, bool is_simple,
                        ExprPtr &expr, unsigned int &nflops,
                        MakeChild &&make_child) {
  using namespace libint2::algebra;
  using libint2::algebra::operator*;

  constexpr int x = 0;
  constexpr int y = 1;
  constexpr int z = 2;

  F Dx_a{a};
  Dx_a.deriv().inc(x);
  F Dx_b{b};
  Dx_b.deriv().inc(x);
  F Dy_a{a};
  Dy_a.deriv().inc(y);
  F Dy_b{b};
  Dy_b.deriv().inc(y);
  F Dz_a{a};
  Dz_a.deriv().inc(z);
  F Dz_b{b};
  Dz_b.deriv().inc(z);

  switch (q) {
    case 0: {
      auto Dx_a_O_Dx_b = make_child(Dx_a, Dx_b);
      auto Dy_a_O_Dy_b = make_child(Dy_a, Dy_b);
      auto Dz_a_O_Dz_b = make_child(Dz_a, Dz_b);
      if (is_simple) {
        expr = Dx_a_O_Dx_b + Dy_a_O_Dy_b + Dz_a_O_Dz_b;
        nflops += 2;
      }
    } break;
    case 1: {
      auto Dy_a_O_Dz_b = make_child(Dy_a, Dz_b);
      auto Dz_a_O_Dy_b = make_child(Dz_a, Dy_b);
      if (is_simple) {
        expr = Dy_a_O_Dz_b - Dz_a_O_Dy_b;
        nflops += 1;
      }
    } break;
    case 2: {
      auto Dz_a_O_Dx_b = make_child(Dz_a, Dx_b);
      auto Dx_a_O_Dz_b = make_child(Dx_a, Dz_b);
      if (is_simple) {
        expr = Dz_a_O_Dx_b - Dx_a_O_Dz_b;
        nflops += 1;
      }
    } break;
    case 3: {
      auto Dx_a_O_Dy_b = make_child(Dx_a, Dy_b);
      auto Dy_a_O_Dx_b = make_child(Dy_a, Dx_b);
      if (is_simple) {
        expr = Dx_a_O_Dy_b - Dy_a_O_Dx_b;
        nflops += 1;
      }
    } break;
    default:
      throw std::runtime_error("sigma_p_1body_fold: invalid quaternionic index");
  }
}

}  // namespace libint2

#endif  // LIBINT_RKB_FOLD_1BODY_H
