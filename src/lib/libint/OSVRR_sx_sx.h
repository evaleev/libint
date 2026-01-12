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

#ifndef _libint2_src_lib_libint_osvrrsxsx_h_
#define _libint2_src_lib_libint_osvrrsxsx_h_

#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <util_types.h>

#include <cstdlib>

#ifdef __GNUC__
#pragma implementation
#endif

namespace libint2 {

template <int part, int Lb, int Ld, bool unit_a, bool vectorize>
struct OSVRR_sx_sx {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4);
};

/** builds (0b|0d)^(m)
    src0 = (0b-1|0d)^(m) // ignored if unit_a = true
    src1 = (0b-1|0d)^(m+1)
    src2 = (0b-2|0d)^(m)
    src3 = (0b-2|0d)^(m+1)
    src4 = (0b-1|0d-1)^(m+1)
 **/
template <int Lb, int Ld, bool unit_a, bool vectorize>
struct OSVRR_sx_sx<0, Lb, Ld, unit_a, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src0,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src2,
                      const LIBINT2_REALTYPE* src3,
                      const LIBINT2_REALTYPE* src4) {
    // works for (sd|sp) and higher
    assert(not(Lb < 2 || Ld < 1));

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    const unsigned int Nd = INT_NCART(Ld);
    const unsigned int NdV = Nd * veclen;

    int bx, by, bz;
    FOR_CART(bx, by, bz, Lb)

    int b[3];
    b[0] = bx;
    b[1] = by;
    b[2] = bz;

    enum XYZ { x = 0, y = 1, z = 2 };
    // Build along x, if possible
    XYZ xyz = z;
    if (by != 0) xyz = y;
    if (bx != 0) xyz = x;
    --b[xyz];

    // redirect
    const LIBINT2_REALTYPE *PB, *WP;
    switch (xyz) {
      case x:
#if LIBINT2_DEFINED(eri, PB_x)
        if (not unit_a) PB = inteval->PB_x;
#endif
        WP = inteval->WP_x;
        break;
      case y:
#if LIBINT2_DEFINED(eri, PB_y)
        if (not unit_a) PB = inteval->PB_y;
#endif
        WP = inteval->WP_y;
        break;
      case z:
#if LIBINT2_DEFINED(eri, PB_z)
        if (not unit_a) PB = inteval->PB_z;
#endif
        WP = inteval->WP_z;
        break;
    }

    const unsigned int ibm1 = INT_CARTINDEX(Lb - 1, b[0], b[1]);
    const unsigned int bm10d0_offset = ibm1 * NdV;
    const LIBINT2_REALTYPE* src0_ptr = unit_a ? 0 : src0 + bm10d0_offset;
    const LIBINT2_REALTYPE* src1_ptr = src1 + bm10d0_offset;

    // if b-2_xyz exists, include (0 b-2_xyz | 0 d)
    if (b[xyz] > 0) {
      --b[xyz];
      const unsigned int ibm2 = INT_CARTINDEX(Lb - 2, b[0], b[1]);
      const unsigned int bm20d0_offset = ibm2 * NdV;
      ++b[xyz];
      const LIBINT2_REALTYPE* src2_ptr = src2 + bm20d0_offset;
      const LIBINT2_REALTYPE* src3_ptr = src3 + bm20d0_offset;
      const LIBINT2_REALTYPE bxyz = (LIBINT2_REALTYPE)b[xyz];

      unsigned int dv = 0;
      for (unsigned int d = 0; d < Nd; ++d) {
        for (unsigned int v = 0; v < veclen; ++v, ++dv) {
          LIBINT2_REALTYPE value =
              WP[v] * src1_ptr[dv] +
              bxyz * inteval->oo2z[v] *
                  (src2_ptr[dv] - inteval->roz[v] * src3_ptr[dv]);
          if (not unit_a) value += PB[v] * src0_ptr[dv];
          target[dv] = value;
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += (unit_a ? 6 : 8) * NdV;
#endif

    } else {
      unsigned int dv = 0;
      for (unsigned int d = 0; d < Nd; ++d) {
        for (unsigned int v = 0; v < veclen; ++v, ++dv) {
          LIBINT2_REALTYPE value = WP[v] * src1_ptr[dv];
          if (not unit_a) value += PB[v] * src0_ptr[dv];
          target[dv] = value;
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += (unit_a ? 1 : 3) * NdV;
#endif
    }

    {
      const unsigned int Ndm1 = INT_NCART(Ld - 1);
      const unsigned int Ndm1V = Ndm1 * veclen;
      const unsigned int bm10dm10_offset = ibm1 * Ndm1V;
      const LIBINT2_REALTYPE* src4_ptr = src4 + bm10dm10_offset;

      // loop over d-1 shell and include (0 b-1_xyz | 0 d-1_xyz) to (0 b | 0 d)
      int dx, dy, dz;
      FOR_CART(dx, dy, dz, Ld - 1)

      int d[3];
      d[0] = dx;
      d[1] = dy;
      d[2] = dz;
      ++d[xyz];

      const unsigned int dc = INT_CARTINDEX(Ld, d[0], d[1]);
      const unsigned int dc_offset = dc * veclen;
      LIBINT2_REALTYPE* tptr = target + dc_offset;
      const LIBINT2_REALTYPE dxyz = (LIBINT2_REALTYPE)d[xyz];
      for (unsigned int v = 0; v < veclen; ++v) {
        tptr[v] += dxyz * inteval->oo2ze[v] * src4_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src4_ptr += veclen;

      END_FOR_CART
    }

    target += NdV;

    END_FOR_CART
  }
};

#if 0
  /** builds (0b|0d)^(m)
      src0 = (0b|0d-1)^(m)
      src1 = (0b|0d-1)^(m+1)
      src2 = (0b|0d-2)^(m)
      src3 = (0b|0d-2)^(m+1)
      src4 = (0b-1|0d-1)^(m+1)
   */
  template <int Lb, int Ld, bool vectorize> struct OSVRR_sx_sx<1,Lb,Ld,vectorize> {

    static void compute(const Libint_t* inteval,
        LIBINT2_REALTYPE* target,
        const LIBINT2_REALTYPE* src0,
        const LIBINT2_REALTYPE* src1,
        const LIBINT2_REALTYPE* src2,
        const LIBINT2_REALTYPE* src3,
        const LIBINT2_REALTYPE* src4) {

      // works for (sp|sd) and higher
      if (Lb < 1 || Ld < 2)
        abort();

      // ALGORITHM
      // loop over functions in d
      //   decide in which direction can build (x, y, or z)
      //   decide whether d-2 exists
      //   loop over functions in b
      //     include contributions from (0 b | 0 d-1)^(m), (0 b | 0 d-1)^(m+1),
      //       and possibly (0 b | 0 d-2)^(m), (0 b | 0 d-2)^(m+1) for each b
      //   end of loop over b
      //   loop over b-1
      //     include contribution from (0 b-1 | 0 d-1)^(m+1)
      //   end of loop over b-1
      // end of loop over d

      const unsigned int veclen = vectorize ? inteval->veclen : 1;

      const unsigned int Nb = INT_NCART(Lb);
      const unsigned int Nd = INT_NCART(Ld);
      const unsigned int Ndv = Nd * veclen;
      const unsigned int Ndm1 = INT_NCART(Ld-1);
      const unsigned int Ndm1v = Ndm1 * veclen;
      const unsigned int Ndm2 = INT_NCART(Ld-2);
      const unsigned int Ndm2v = Ndm2 * veclen;

      int dx, dy, dz;
      int id = 0;
      FOR_CART(dx, dy, dz, Ld)

        int d[3]; d[0] = dx;  d[1] = dy;  d[2] = dz;

        enum XYZ {x=0, y=1, z=2};
        // Build along x, if possible
        XYZ xyz = z;
        if (dy != 0) xyz = y;
        if (dx != 0) xyz = x;
        --d[xyz];

        // redirect
        const LIBINT2_REALTYPE *QD, *WQ;
        switch(xyz) {
          case x:
            QD = inteval->QD_x;
            WQ = inteval->WQ_x;
            break;
          case y:
            QD = inteval->QD_y;
            WQ = inteval->WQ_y;
            break;
          case z:
            QD = inteval->QD_z;
            WQ = inteval->WQ_z;
            break;
        }

        const unsigned int idm1 = INT_CARTINDEX(Ld-1,d[0],d[1]);
        const unsigned int d0_offset = id * veclen;
        const unsigned int dm10_offset = idm1 * veclen;
        LIBINT2_REALTYPE* target_ptr = target + d0_offset;
        const LIBINT2_REALTYPE* src0_ptr = src0 + dm10_offset;
        const LIBINT2_REALTYPE* src1_ptr = src1 + dm10_offset;

        // if d-2_xyz exists, include (0 b | 0 d-2_xyz)
        if (d[xyz] > 0) {
          --d[xyz];
          const unsigned int idm2 = INT_CARTINDEX(Ld-2,d[0],d[1]);
          const unsigned int dm20_offset = idm2 * veclen;
          ++d[xyz];
          const LIBINT2_REALTYPE* src2_ptr = src2 + dm20_offset;
          const LIBINT2_REALTYPE* src3_ptr = src3 + dm20_offset;
          const LIBINT2_REALTYPE dxyz = (LIBINT2_REALTYPE)d[xyz];

          for(unsigned int b = 0; b < Nb; ++b) {
            for(unsigned int v=0; v<veclen; ++v) {
              target_ptr[v] = QD[v] * src0_ptr[v] + WQ[v] * src1_ptr[v]
                            + dxyz * inteval->oo2e[v] * (src2_ptr[v] - inteval->roe[v] * src3_ptr[v]);
            }
            target_ptr += Ndv;
            src0_ptr += Ndm1v;
            src1_ptr += Ndm1v;
            src2_ptr += Ndm2v;
            src3_ptr += Ndm2v;
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 8 * Nb * veclen;
#endif

        }
        else {
          for(unsigned int b = 0; b < Nb; ++b) {
            for(unsigned int v=0; v<veclen; ++v) {
              target_ptr[v] = QD[v] * src0_ptr[v] + WQ[v] * src1_ptr[v];
            }
            target_ptr += Ndv;
            src0_ptr += Ndm1v;
            src1_ptr += Ndm1v;
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 3 * Nb * veclen;
#endif
        }

        {
          const LIBINT2_REALTYPE* src4_ptr = src4 + dm10_offset;

          // loop over b-1 shell and include (0 b-1_xyz | 0 d-1_xyz) to (0 b | 0 d)
          int bx, by, bz;
          FOR_CART(bx, by, bz, Lb-1)

            int b[3]; b[0] = bx;  b[1] = by;  b[2] = bz;
            ++b[xyz];

            const unsigned int ib = INT_CARTINDEX(Lb,b[0],b[1]);
            const unsigned int b0d0_offset = ib * Ndv + d0_offset;
            LIBINT2_REALTYPE* target_ptr = target + b0d0_offset;
            const LIBINT2_REALTYPE bxyz = (LIBINT2_REALTYPE)b[xyz];
            for(unsigned int v=0; v<veclen; ++v) {
              target_ptr[v] += bxyz * inteval->oo2ze[v] * src4_ptr[v];
            }
#if LIBINT2_FLOP_COUNT
            inteval->nflops[0] += 3 * veclen;
#endif
            src4_ptr += Ndm1v;

          END_FOR_CART // end of loop over b
        }

        ++id;

      END_FOR_CART // end of loop over d

    }

  };
#endif

template <int part, int Lb, int Ld, bool vectorize>
struct OSAVRR_sx_sx {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src4);
};

/** builds (0b|0d)^(m)
    src1 = (0b-1|0d)^(m+1)
    src4 = (0b-1|0d-1)^(m+1)
 **/
template <int Lb, int Ld, bool vectorize>
struct OSAVRR_sx_sx<0, Lb, Ld, vectorize> {
  static void compute(const Libint_t* inteval, LIBINT2_REALTYPE* target,
                      const LIBINT2_REALTYPE* src1,
                      const LIBINT2_REALTYPE* src4) {
    // works for (sd|sp) and higher
    assert(not(Lb < 2 || Ld < 1));

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    const unsigned int Nd = INT_NCART(Ld);
    const unsigned int NdV = Nd * veclen;

    int bx, by, bz;
    FOR_CART(bx, by, bz, Lb)

    int b[3];
    b[0] = bx;
    b[1] = by;
    b[2] = bz;

    enum XYZ { x = 0, y = 1, z = 2 };
    // Build along x, if possible
    XYZ xyz = z;
    if (by != 0) xyz = y;
    if (bx != 0) xyz = x;
    --b[xyz];

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

    const unsigned int ibm1 = INT_CARTINDEX(Lb - 1, b[0], b[1]);
    const unsigned int bm10d0_offset = ibm1 * NdV;
    const LIBINT2_REALTYPE* src1_ptr = src1 + bm10d0_offset;

    {
      unsigned int dv = 0;
      for (unsigned int d = 0; d < Nd; ++d) {
        for (unsigned int v = 0; v < veclen; ++v, ++dv) {
          target[dv] = WP[v] * src1_ptr[dv];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += NdV;
#endif
    }

    {
      const unsigned int Ndm1 = INT_NCART(Ld - 1);
      const unsigned int Ndm1V = Ndm1 * veclen;
      const unsigned int bm10dm10_offset = ibm1 * Ndm1V;
      const LIBINT2_REALTYPE* src4_ptr = src4 + bm10dm10_offset;

      // loop over d-1 shell and include (0 b-1_xyz | 0 d-1_xyz) to (0 b | 0 d)
      int dx, dy, dz;
      FOR_CART(dx, dy, dz, Ld - 1)

      int d[3];
      d[0] = dx;
      d[1] = dy;
      d[2] = dz;
      ++d[xyz];

      const unsigned int dc = INT_CARTINDEX(Ld, d[0], d[1]);
      const unsigned int dc_offset = dc * veclen;
      LIBINT2_REALTYPE* tptr = target + dc_offset;
      const LIBINT2_REALTYPE dxyz = (LIBINT2_REALTYPE)d[xyz];
      for (unsigned int v = 0; v < veclen; ++v) {
        tptr[v] += dxyz * inteval->oo2ze[v] * src4_ptr[v];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 3 * veclen;
#endif
      src4_ptr += veclen;

      END_FOR_CART
    }

    target += NdV;

    END_FOR_CART
  }
};

};  // namespace libint2

#endif  // header guard
