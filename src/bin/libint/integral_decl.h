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

#ifndef _libint2_src_bin_libint_integraldecl_h_
#define _libint2_src_bin_libint_integraldecl_h_

namespace libint2 {

template <class Oper, class BFS, class BraSetType, class KetSetType,
          class AuxQuanta>
class GenIntegralSet;
#if LIBINT_SUPPORT_ONEBODYINTS
template <class Oper, class BFS, class AuxQuanta>
class GenIntegralSet_1_1;
#endif  // LIBINT_SUPPORT_ONEBODYINTS
template <class Oper, class BFS, class AuxQuanta>
class GenIntegralSet_11_11;

#if 0
  template <class BFS> class TwoPRep_11_11;
#endif
template <class BFS, int K>
class R12kG12_11_11;
template <class BFS, int K>
class TiG12_11_11;
template <class BFS>
class R1dotR1G12_11_11;
template <class BFS>
class R2dotR2G12_11_11;
template <class BFS>
class R1dotR2G12_11_11;

};  // namespace libint2

#endif
