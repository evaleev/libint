/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_include_braket_h_
#define _libint2_include_braket_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "The simple Libint API requires C++11 support"
#endif

namespace libint2 {
/// types of shell sets supported by Engine, in chemist notation (i.e. '_'
/// separates particles)
/// \warning macro \c BOOST_PP_NBODY_BRAKET_RANK_TUPLE include the ranks of all
/// brakets in \c BraKet
///          and macro \c BOOST_PP_NBODY_BRAKET_MAX_INDEX must be equal to the
///          max value in this enum
enum class BraKet {
  x_x = 0,
  xx_xx,
  xs_xx,
  xx_xs,
  xs_xs,
  invalid = -1,
  first_1body_braket = x_x,
  last_1body_braket = x_x,
  first_2body_braket = xx_xx,
  last_2body_braket = xs_xs,
  first_braket = first_1body_braket,
  last_braket = last_2body_braket
};

}  // namespace libint2

#endif /* _libint2_include_braket_h_*/
