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

#ifndef _libint2_src_lib_libint_itrxsxs_h_
#define _libint2_src_lib_libint_itrxsxs_h_

#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <util_types.h>

#include <cstdlib>

namespace libint2 {

template <int part, int La, int Lc, bool InBra, bool vectorize>
struct ITR_xs_xs {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3);
};

/** builds (a 0|c0) from
    src0 = (a-1 0|c 0)
    src1 = (a-1 0|c+1 0)
    src2 = (a-2 0|c 0)
    src3 = (a-1 0|c-1 0)
 **/
template <int La, int Lc, bool InBra, bool vectorize>
struct ITR_xs_xs<0, La, Lc, InBra, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3) {
    // works for (ds|ps) and higher
    if (La < 2 || Lc < 1) abort();

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
    const LIBINT2_REALTYPE* pfac0;
    switch (xyz) {
      case x:
        pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_x)
                      inteval->TwoPRepITR_pfac0_0_0_x
#else
                      0
#endif
                      :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
                      inteval->TwoPRepITR_pfac0_0_1_x;
#else
                      0
#endif
        ;
        break;
      case y:
        pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_y)
                      inteval->TwoPRepITR_pfac0_0_0_y
#else
                      0
#endif
                      :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_y)
                      inteval->TwoPRepITR_pfac0_0_1_y;
#else
                      0
#endif
        ;
        break;
      case z:
        pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_z)
                      inteval->TwoPRepITR_pfac0_0_0_z
#else
                      0
#endif
                      :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_z)
                      inteval->TwoPRepITR_pfac0_0_1_z;
#else
                      0
#endif
        ;
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
      const LIBINT2_REALTYPE* src2_ptr = src2 + am20c0_offset;
      const LIBINT2_REALTYPE axyz = (LIBINT2_REALTYPE)a[xyz];

      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] =
              pfac0[v] * src0_ptr[cv] + axyz * inteval->oo2z[v] * src2_ptr[cv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 4 * NcV;
#endif

    } else {
      unsigned int cv = 0;
      for (unsigned int c = 0; c < Nc; ++c) {
        for (unsigned int v = 0; v < veclen; ++v, ++cv) {
          target[cv] = pfac0[v] * src0_ptr[cv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 1 * NcV;
#endif
    }

    {
      const unsigned int Ncp1 = INT_NCART(Lc + 1);
      const unsigned int Ncp1V = Ncp1 * veclen;
      const unsigned int am10cp10_offset = iam1 * Ncp1V;
      const LIBINT2_REALTYPE* src1_ptr = src1 + am10cp10_offset;

      // loop over c shell and include (a-1_xyz 0 | c+1_xyz 0) to (a 0 | c 0)
      int cx, cy, cz;
      int cv = 0;
      FOR_CART(cx, cy, cz, Lc)

      int c[3];
      c[0] = cx;
      c[1] = cy;
      c[2] = cz;
      ++c[xyz];

      const unsigned int cp1 = INT_CARTINDEX(Lc + 1, c[0], c[1]);
      const unsigned int cp1_offset = cp1 * veclen;
      const LIBINT2_REALTYPE* sptr = src1_ptr + cp1_offset;
      const LIBINT2_REALTYPE cxyz = (LIBINT2_REALTYPE)c[xyz];
      for (unsigned int v = 0; v < veclen; ++v, ++cv) {
        target[cv] -= inteval->eoz[v] * sptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 2 * veclen;
#endif

      END_FOR_CART
    }

    {
      const unsigned int Ncm1 = INT_NCART(Lc - 1);
      const unsigned int Ncm1V = Ncm1 * veclen;
      const unsigned int am10cm10_offset = iam1 * Ncm1V;
      const LIBINT2_REALTYPE* src3_ptr = src3 + am10cm10_offset;

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
        tptr[v] += cxyz * inteval->oo2z[v] * src3_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src3_ptr += veclen;

      END_FOR_CART
    }

    target += NcV;

    END_FOR_CART
  }
};

#if 0
  /** builds (a 0|c0) from
      src0 = (a 0|c-1 0)
      src1 = (a+1 0|c-1 0)
      src2 = (a 0|c-2 0)
      src3 = (a-1 0|c-1 0)
   **/
  template <int La, int Lc, bool InBra, bool vectorize> struct ITR_xs_xs<1,La,Lc,InBra,vectorize> {

    static void compute(const Libint_t* inteval,
        LIBINT2_REALTYPE* target,
        const LIBINT2_REALTYPE* src0,
        const LIBINT2_REALTYPE* src1,
        const LIBINT2_REALTYPE* src2,
        const LIBINT2_REALTYPE* src3) {

      // works for (ps|ds) and higher
      if (La < 1 || Lc < 2)
        abort();

      const unsigned int veclen = vectorize ? inteval->veclen : 1;

      const unsigned int Na = INT_NCART(La);
      const unsigned int NaV = Na * veclen;

      int cx, cy, cz;
      FOR_CART(cx, cy, cz, Lc)

        int c[3]; c[0] = cx;  c[1] = cy;  c[2] = cz;

        enum XYZ {x=0, y=1, z=2};
        // Build along x, if possible
        XYZ xyz = z;
        if (cy != 0) xyz = y;
        if (cx != 0) xyz = x;
        --c[xyz];

        // redirect
        const LIBINT2_REALTYPE *pfac0;
        switch(xyz) {
          case x:
            pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_x)
                inteval->TwoPRepITR_pfac0_1_0_x
#else
                0
#endif
                :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
                inteval->TwoPRepITR_pfac0_0_1_x;
#else
                0
#endif
                ;
            break;
          case y:
            pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_y)
                inteval->TwoPRepITR_pfac0_1_0_y
#else
                0
#endif
                :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_y)
                inteval->TwoPRepITR_pfac0_1_1_y;
#else
                0
#endif
                ;
            break;
          case z:
            pfac0 = InBra ?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_z)
                inteval->TwoPRepITR_pfac0_1_0_z
#else
                0
#endif
                :
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_z)
                inteval->TwoPRepITR_pfac0_1_1_z;
#else
                0
#endif
                ;
            break;
        }

        const unsigned int icm1 = INT_CARTINDEX(Lc-1,c[0],c[1]);
        const unsigned int a0cm10_offset = icm1 * NaV;
        const LIBINT2_REALTYPE* src0_ptr = src0 + a0cm10_offset;

        // if c-2_xyz exists, include (a 0 | c-2_xyz 0)
        if (c[xyz] > 0) {
          --c[xyz];
          const unsigned int icm2 = INT_CARTINDEX(Lc-2,c[0],c[1]);
          const unsigned int a0cm20_offset = icm2 * NaV;
          ++c[xyz];
          const LIBINT2_REALTYPE* src2_ptr = src2 + am20c0_offset;
          const LIBINT2_REALTYPE axyz = (LIBINT2_REALTYPE)c[xyz];

          unsigned int cv = 0;
          for(unsigned int c = 0; c < Nc; ++c) {
            for(unsigned int v=0; v<veclen; ++v, ++cv) {
              target[cv] = pfac0[v] * src0_ptr[cv] + axyz * inteval->oo2z[v] * src2_ptr[cv];
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 4 * NcV;
#endif

        }
        else {
          unsigned int cv = 0;
          for(unsigned int c = 0; c < Nc; ++c) {
            for(unsigned int v=0; v<veclen; ++v, ++cv) {
              target[cv] = pfac0[v] * src0_ptr[cv];
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 1 * NcV;
#endif
        }

        {
          const unsigned int Ncp1 = INT_NCART(Lc+1);
          const unsigned int Ncp1V = Ncp1 * veclen;
          const unsigned int am10cp10_offset = iam1 * Ncp1V;
          const LIBINT2_REALTYPE* src1_ptr = src1 + am10cp10_offset;

          // loop over c shell and include (c-1_xyz 0 | c+1_xyz 0) to (c 0 | c 0)
          int cx, cy, cz;
          int cv = 0;
          FOR_CART(cx, cy, cz, Lc)

            int c[3]; c[0] = cx;  c[1] = cy;  c[2] = cz;
            ++c[xyz];

            const unsigned int cp1 = INT_CARTINDEX(Lc+1,c[0],c[1]);
            const unsigned int cp1_offset = cp1 * veclen;
            const LIBINT2_REALTYPE* sptr = src1_ptr + cp1_offset;
            const LIBINT2_REALTYPE cxyz = (LIBINT2_REALTYPE)c[xyz];
            for(unsigned int v=0; v<veclen; ++v, ++cv) {
              target[cv] -= inteval->eoz[v] * sptr[v];
            }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 2 * veclen;
#endif

          END_FOR_CART
        }


        {
          const unsigned int Ncm1 = INT_NCART(Lc-1);
          const unsigned int Ncm1V = Ncm1 * veclen;
          const unsigned int am10cm10_offset = iam1 * Ncm1V;
          const LIBINT2_REALTYPE* src3_ptr = src3 + am10cm10_offset;

          // loop over c-1 shell and include (c-1_xyz 0 | c-1_xyz 0) to (c 0 | c 0)
          int cx, cy, cz;
          FOR_CART(cx, cy, cz, Lc-1)

            int c[3]; c[0] = cx;  c[1] = cy;  c[2] = cz;
            ++c[xyz];

            const unsigned int cc = INT_CARTINDEX(Lc,c[0],c[1]);
            const unsigned int cc_offset = cc * veclen;
            LIBINT2_REALTYPE* tptr = target + cc_offset;
            const LIBINT2_REALTYPE cxyz = (LIBINT2_REALTYPE)c[xyz];
            for(unsigned int v=0; v<veclen; ++v) {
              tptr[v] += cxyz * inteval->oo2z[v] * src3_ptr[v];
            }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 3 * veclen;
#endif
            src3_ptr += veclen;

          END_FOR_CART
        }

        target += NcV;

      END_FOR_CART

    }

  };
#endif

};  // namespace libint2

#endif  // header guard
