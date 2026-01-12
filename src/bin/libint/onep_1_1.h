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

#ifndef _libint2_src_bin_libint_onep11_h_
#define _libint2_src_bin_libint_onep11_h_

#include <integral.h>
#include <integral_1_1.h>

namespace libint2 {

/// (s|V|s) (V=electrostatic potential operator) shell quartet is not
/// precomputed, but the integral is
template <>
inline bool GenIntegralSet_1_1<CGF, ElecPotOper, mType>::this_precomputed()
    const {
  /// uncontracted (s|s) are precomputed
  if (parent_type::bra_.member(0, 0).zero() &&
      parent_type::ket_.member(0, 0).zero() &&
      !parent_type::bra_.member(0, 0).contracted() &&
      !parent_type::ket_.member(0, 0).contracted() &&
      parent_type::bra_.member(0, 0).deriv().zero() &&
      parent_type::ket_.member(0, 0).deriv().zero())
    return true;
  else
    return false;
}

/// always unroll (s|V|s)^(m)
template <>
inline bool GenIntegralSet_1_1<CGShell, ElecPotOper, mType>::auto_unroll()
    const {
  /// uncontracted (s|s) are precomputed
  if (parent_type::bra_.member(0, 0).zero() &&
      parent_type::ket_.member(0, 0).zero() &&
      !parent_type::bra_.member(0, 0).contracted() &&
      !parent_type::ket_.member(0, 0).contracted() &&
      parent_type::bra_.member(0, 0).deriv().zero() &&
      parent_type::ket_.member(0, 0).deriv().zero())
    return true;
  else
    return false;
}

};  // namespace libint2

#endif
