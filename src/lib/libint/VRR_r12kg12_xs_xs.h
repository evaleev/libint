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

#ifndef _libint2_src_lib_libint_vrrr12kg12xsxs_h_
#define _libint2_src_lib_libint_vrrr12kg12xsxs_h_

#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <util_types.h>

#include <cstdlib>

namespace libint2 {

template <int part, int La, int Lc, int K, bool vectorize>
struct VRR_r12kg12_xs_xs {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4,
                      const LIBINT2_REALTYPE* src5);
};

/** builds (a0| G_K |c0), where G_K = r12^K * G12, for K >= 0. For K == -1 use
 OSVRR (just make sure prefactors are modified accordingly). src0 =
 (a-10|G_K|c0) src1 = (a-20|G_K|c0) src2 = (a-10|G_K|c-10) src3 = (a0|G_K-2|c0)
    src4 = (a-10|G_K-2|c+10)
    src5 = (a-10|G_K-2|c0)

    this code will also build (0a| G_K |0c) as long as the RR prefactors are
 defined appropriately.
 **/
template <int La, int Lc, int K, bool vectorize>
struct VRR_r12kg12_xs_xs<0, La, Lc, K, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4,
                      const LIBINT2_REALTYPE* src5) {
    // works for (ds|ps) and higher, K >= 0
    if ((La < 2 && Lc < 1) || K < 0) abort();

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    const unsigned int Nc = INT_NCART(Lc);
    const unsigned int NcV = Nc * veclen;
    const unsigned int Ncm1 = INT_NCART(Lc - 1);
    const unsigned int Ncm1V = Ncm1 * veclen;
    const unsigned int Ncp1 = INT_NCART(Lc + 1);
    const unsigned int Ncp1V = Ncp1 * veclen;

    int ax, ay, az;
    FOR_CART(ax, ay, az, La)

    int a[3];
    a[0] = ax;
    a[1] = ay;
    a[2] = az;

    enum XYZ { x = 0, y = 1, z = 2 };
    // Build along x, if possible
    XYZ xyz = z;
    if (ay != 0) xyz = y;
    if (ax != 0) xyz = x;
    --a[xyz];

    // redirect
    const LIBINT2_REALTYPE *pfac0_0, *pfac4_0;
    switch (xyz) {
      case x:
        pfac0_0 = inteval->R12kG12_pfac0_0_x;
        pfac4_0 = inteval->R12kG12_pfac4_0_x;
        break;
      case y:
        pfac0_0 = inteval->R12kG12_pfac0_0_y;
        pfac4_0 = inteval->R12kG12_pfac4_0_y;
        break;
      case z:
        pfac0_0 = inteval->R12kG12_pfac0_0_z;
        pfac4_0 = inteval->R12kG12_pfac4_0_z;
        break;
    }

    const unsigned int iam1 = INT_CARTINDEX(La - 1, a[0], a[1]);
    const unsigned int am10c0_offset = iam1 * NcV;
    const LIBINT2_REALTYPE* src0_ptr = src0 + am10c0_offset;

    // if a-2_xyz exists, include (a-2_xyz 0 | c 0)
    if (a[xyz] > 0) {
      --a[xyz];
      const unsigned int iam2 = INT_CARTINDEX(La - 2, a[0], a[1]);
      const unsigned int am20c0_offset = iam2 * NcV;
      ++a[xyz];
      const LIBINT2_REALTYPE* src1_ptr = src1 + am20c0_offset;
      const LIBINT2_REALTYPE axyz = (LIBINT2_REALTYPE)a[xyz];

      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] = pfac0_0[v] * src0_ptr[cv] +
                       axyz * inteval->R12kG12_pfac1_0[v] * src1_ptr[cv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 4 * NcV;
#endif

    } else {
      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] = pfac0_0[v] * src0_ptr[cv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += NcV;
#endif
    }

    {
      const unsigned int am10cm10_offset = iam1 * Ncm1V;
      const LIBINT2_REALTYPE* src2_ptr = src2 + am10cm10_offset;

      // loop over c-1 shell and include (a-1_xyz 0 | c-1_xyz 0) to (a 0 | c 0)
      int cx, cy, cz;
      FOR_CART(cx, cy, cz, Lc - 1)

      int c[3];
      c[0] = cx;
      c[1] = cy;
      c[2] = cz;
      ++c[xyz];

      const unsigned int cc = INT_CARTINDEX(Lc, c[0], c[1]);
      const unsigned int cc_offset = cc * veclen;
      LIBINT2_REALTYPE* tptr = target + cc_offset;
      const LIBINT2_REALTYPE cxyz = (LIBINT2_REALTYPE)c[xyz];
      for (unsigned int v = 0; v < veclen; ++v) {
        tptr[v] += cxyz * inteval->R12kG12_pfac2[v] * src2_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src2_ptr += veclen;

      END_FOR_CART
    }

    if (K > 0) {
      ++a[xyz];
      const unsigned int ia = INT_CARTINDEX(La, a[0], a[1]);
      const unsigned int a0c0_offset = ia * NcV;
      --a[xyz];
      const LIBINT2_REALTYPE* src3_ptr = src3 + a0c0_offset;
      const LIBINT2_REALTYPE* src5_ptr = src5 + am10c0_offset;
      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] += K * inteval->R12kG12_pfac3_0[v] *
                        (src3_ptr[cv] + pfac4_0[v] * src5_ptr[cv]);
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 5 * NcV;
#endif

      // loop over c shell and include (a-1_xyz 0 | c+1_xyz 0) to (a 0 | c 0)
      int cx, cy, cz;
      cv = 0;
      const unsigned int am10cp10_offset = iam1 * Ncp1V;
      const LIBINT2_REALTYPE* src4_ptr0 = src4 + am10cp10_offset;
      FOR_CART(cx, cy, cz, Lc)

      int c[3];
      c[0] = cx;
      c[1] = cy;
      c[2] = cz;
      ++c[xyz];

      const unsigned int cc = INT_CARTINDEX(Lc + 1, c[0], c[1]);
      const unsigned int cc_offset = cc * veclen;
      const LIBINT2_REALTYPE* src4_ptr = src4_ptr0 + cc_offset;
      for (unsigned int v = 0; v < veclen; ++v, ++cv) {
        target[cv] -= K * inteval->R12kG12_pfac3_0[v] * src4_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif

      END_FOR_CART
    }

    target += NcV;

    END_FOR_CART

    /** Number of flops = ??? */
    // inteval->nflops[0] = inteval->nflops[0] + 222 * 1 * 1 * veclen;
  }
};

};  // namespace libint2

#endif  // header guard
