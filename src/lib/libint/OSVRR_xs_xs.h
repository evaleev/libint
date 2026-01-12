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

#ifndef _libint2_src_lib_libint_osvrrxsxs_h_
#define _libint2_src_lib_libint_osvrrxsxs_h_

#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <util_types.h>

#include <cstdlib>

namespace libint2 {

template <int part, int La, int Lc, bool unit_b, bool vectorize>
struct OSVRR_xs_xs {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4);
};

/** builds (a 0|c0)^(m)
    src0 = (a-10|c0)^(m) // ignored if unit_b is true
    src1 = (a-10|c0)^(m+1)
    src2 = (a-20|c0)^(m)
    src3 = (a-20|c0)^(m+1)
    src4 = (a-10|c-10)^(m+1)
 **/
template <int La, int Lc, bool unit_b, bool vectorize>
struct OSVRR_xs_xs<0, La, Lc, unit_b, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4) {
    // works for (ds|ps) and higher
    assert(not(La < 2 || Lc < 1));

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    const unsigned int Nc = INT_NCART(Lc);
    const unsigned int NcV = Nc * veclen;

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
    const LIBINT2_REALTYPE *PA, *WP;
    switch (xyz) {
      case x:
#if LIBINT2_DEFINED(eri, PA_x)
        if (not unit_b) PA = inteval->PA_x;
#endif
        WP = inteval->WP_x;
        break;
      case y:
#if LIBINT2_DEFINED(eri, PA_y)
        if (not unit_b) PA = inteval->PA_y;
#endif
        WP = inteval->WP_y;
        break;
      case z:
#if LIBINT2_DEFINED(eri, PA_z)
        if (not unit_b) PA = inteval->PA_z;
#endif
        WP = inteval->WP_z;
        break;
    }

    const unsigned int iam1 = INT_CARTINDEX(La - 1, a[0], a[1]);
    const unsigned int am10c0_offset = iam1 * NcV;
    const LIBINT2_REALTYPE* src0_ptr = unit_b ? 0 : src0 + am10c0_offset;
    const LIBINT2_REALTYPE* src1_ptr = src1 + am10c0_offset;

    // if a-2_xyz exists, include (a-2_xyz 0 | c 0)
    if (a[xyz] > 0) {
      --a[xyz];
      const unsigned int iam2 = INT_CARTINDEX(La - 2, a[0], a[1]);
      const unsigned int am20c0_offset = iam2 * NcV;
      ++a[xyz];
      const LIBINT2_REALTYPE* src2_ptr = src2 + am20c0_offset;
      const LIBINT2_REALTYPE* src3_ptr = src3 + am20c0_offset;
      const LIBINT2_REALTYPE axyz = (LIBINT2_REALTYPE)a[xyz];

      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          LIBINT2_REALTYPE value =
              WP[v] * src1_ptr[cv] +
              axyz * inteval->oo2z[v] *
                  (src2_ptr[cv] - inteval->roz[v] * src3_ptr[cv]);
          if (not unit_b) value += PA[v] * src0_ptr[cv];
          target[cv] = value;
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += (unit_b ? 6 : 8) * NcV;
#endif

    } else {
      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          LIBINT2_REALTYPE value = WP[v] * src1_ptr[cv];
          if (not unit_b) value += PA[v] * src0_ptr[cv];
          target[cv] = value;
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += (unit_b ? 1 : 3) * NcV;
#endif
    }

    {
      const unsigned int Ncm1 = INT_NCART(Lc - 1);
      const unsigned int Ncm1V = Ncm1 * veclen;
      const unsigned int am10cm10_offset = iam1 * Ncm1V;
      const LIBINT2_REALTYPE* src4_ptr = src4 + am10cm10_offset;

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
        tptr[v] += cxyz * inteval->oo2ze[v] * src4_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src4_ptr += veclen;

      END_FOR_CART
    }

    target += NcV;

    END_FOR_CART

    /** Number of flops = ??? */
    // inteval->nflops[0] = inteval->nflops[0] + 222 * 1 * 1 * veclen;
  }
};

// Ahlrichs' extension of OS VRR
template <int part, int La, int Lc, bool vectorize>
struct OSAVRR_xs_xs {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src4);
};

/** builds (a 0|c0)^(m)
    src1 = (a-10|c0)^(m+1)
    src4 = (a-10|c-10)^(m+1)
 **/
template <int La, int Lc, bool vectorize>
struct OSAVRR_xs_xs<0, La, Lc, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src4) {
    // works for (ps|ps) and higher
    assert(not(La < 1 || Lc < 1));

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    const unsigned int Nc = INT_NCART(Lc);
    const unsigned int NcV = Nc * veclen;

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
    const LIBINT2_REALTYPE* WP;
    switch (xyz) {
      case x:
        WP = inteval->WP_x;
        break;
      case y:
        WP = inteval->WP_y;
        break;
      case z:
        WP = inteval->WP_z;
        break;
    }

    const unsigned int iam1 = INT_CARTINDEX(La - 1, a[0], a[1]);
    const unsigned int am10c0_offset = iam1 * NcV;
    const LIBINT2_REALTYPE* src1_ptr = src1 + am10c0_offset;

    {
      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] = WP[v] * src1_ptr[cv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += NcV;
#endif
    }

    {
      const unsigned int Ncm1 = INT_NCART(Lc - 1);
      const unsigned int Ncm1V = Ncm1 * veclen;
      const unsigned int am10cm10_offset = iam1 * Ncm1V;
      const LIBINT2_REALTYPE* src4_ptr = src4 + am10cm10_offset;

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
        tptr[v] += cxyz * inteval->oo2ze[v] * src4_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src4_ptr += veclen;

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
