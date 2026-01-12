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

#ifndef _libint2_src_lib_libint_genericgaussderiv_h_
#define _libint2_src_lib_libint_genericgaussderiv_h_

#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <util_types.h>

#include <cassert>
#include <cstdlib>

namespace libint2 {

/** builds ( ... d a / d r_dir ... )
    src0 = ( ... a+1 ... )
    src1 = ( ... a-1 ... )
 **/
template <int L,  //!< the angular momentum of this shell (i.e. the shell being
                  //!< differentiated)
          bool vectorize>
struct GenericGaussDeriv {
  static void compute(
      const Libint_t* inteval, LIBINT2_REALTYPE* target,
      const LIBINT2_REALTYPE* src0, const LIBINT2_REALTYPE* src1,
      unsigned int
          highdim,  //!< number of functions more significant than this shell
      unsigned int
          lowdim,  //!<  number of functions less significant than this shell
      unsigned int dir,  //!< Cartesian direction of the derivative (0, 1, 2)
      const LIBINT2_REALTYPE (
          &two_alpha)[LIBINT2_MAX_VECLEN]  //!< The gaussian exponent of the
                                           //!< function being differentiated)
  );
};

};  // namespace libint2

#endif  // header guard
