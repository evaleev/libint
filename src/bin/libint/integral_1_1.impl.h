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

#ifndef _libint2_src_bin_libint_integral11impl_h_
#define _libint2_src_bin_libint_integral11impl_h_

namespace libint2 {

template <class BFS, class Oper, class AuxQuanta>
bool GenIntegralSet_1_1<BFS, Oper, AuxQuanta>::this_precomputed() const {
  return false;
}

template <class BFS, class Oper, class AuxQuanta>
bool GenIntegralSet_1_1<BFS, Oper, AuxQuanta>::auto_unroll() const {
  return false;
}

};  // namespace libint2

#endif
