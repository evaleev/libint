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

#ifndef _libint2_src_lib_libint_genericgaussderivimpl_h_
#define _libint2_src_lib_libint_genericgaussderivimpl_h_

#include <GenericGaussDeriv.h>

namespace libint2 {

/** builds ( ... d a / d r_dir ... )
    src0 = ( ... a+1 ... )
    src1 = ( ... a-1 ... )
 **/
template <int L,  //!< the angular momentum of this shell (i.e. the shell being
                  //!< differentiated)
          bool vectorize>
void GenericGaussDeriv<L, vectorize>::compute(
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
) {
  const unsigned int veclen = vectorize ? inteval->veclen : 1;
  const unsigned int lveclen = lowdim * veclen;

  const unsigned int N = INT_NCART(L);
  const unsigned int Np1 = INT_NCART(L + 1);
  const unsigned int Nm1 =
      L > 0 ? INT_NCART(L - 1) : 0;  // L == 0 case will not use Nm1

  for (unsigned int h = 0, target_idx = 0; h < highdim; ++h) {
    const unsigned int hNp1 = h * Np1;
    const unsigned int hNm1 = h * Nm1;

    int ax, ay, az;
    FOR_CART(ax, ay, az, L)

    int a[3];
    a[0] = ax;
    a[1] = ay;
    a[2] = az;

    // d |a) / dRi = 2 alpha a+1_i
    ++a[dir];
    const unsigned int iap1 = INT_CARTINDEX(L + 1, a[0], a[1]);
    const unsigned int ap1_offset = (hNp1 + iap1) * lveclen;
    const LIBINT2_REALTYPE* src0_ptr = src0 + ap1_offset;
    --a[dir];

    const bool have_am1 = (a[dir] > 0);
    // d |a) / dRi -=  a_i a-1_i
    unsigned int iam1;
    unsigned int am1_offset;
    const LIBINT2_REALTYPE* src1_ptr;
    if (have_am1) {
      --a[dir];
      iam1 = INT_CARTINDEX(L - 1, a[0], a[1]);
      am1_offset = (hNm1 + iam1) * lveclen;
      src1_ptr = src1 + am1_offset;
      ++a[dir];
    }
    LIBINT2_REALTYPE adir_real = LIBINT2_REALTYPE(a[dir]);

    if (have_am1)
      for (unsigned int l = 0, lv = 0; l < lowdim; ++l) {
        for (unsigned int v = 0; v < veclen; ++v, ++lv, ++target_idx) {
          target[target_idx] =
              two_alpha[v] * src0_ptr[lv] - adir_real * src1_ptr[lv];
        }
      }
    else
      for (unsigned int l = 0, lv = 0; l < lowdim; ++l) {
        for (unsigned int v = 0; v < veclen; ++v, ++lv, ++target_idx) {
          target[target_idx] = two_alpha[v] * src0_ptr[lv];
        }
      }

#if LIBINT2_FLOP_COUNT
    inteval->nflops[0] += (have_am1 ? 3 : 1) * lveclen;
#endif

    END_FOR_CART  // end of loop over a
  }
}

};  // namespace libint2

#endif  // header guard
