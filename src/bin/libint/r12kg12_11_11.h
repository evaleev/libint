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

#ifndef _libint2_src_bin_libint_r12kg121111_h_
#define _libint2_src_bin_libint_r12kg121111_h_

#include <integral.h>
#include <integral_11_11.h>

namespace libint2 {

template <>
inline bool GenIntegralSet_11_11<CGF, R12kG12, mType>::this_precomputed()
    const {
  if (parent_type::bra_.member(0, 0).zero() &&
      parent_type::bra_.member(1, 0).zero() &&
      parent_type::ket_.member(0, 0).zero() &&
      parent_type::ket_.member(1, 0).zero())
    return true;
  else
    return false;
}

template <>
inline bool GenIntegralSet_11_11<CGShell, R12kG12, mType>::auto_unroll() const {
  if (parent_type::bra_.member(0, 0).zero() &&
      parent_type::bra_.member(1, 0).zero() &&
      parent_type::ket_.member(0, 0).zero() &&
      parent_type::ket_.member(1, 0).zero())
    return true;
  else
    return false;
}
};  // namespace libint2

#endif
