
#ifndef _libint2_src_lib_libint_osvrrsxsx_h_
#define _libint2_src_lib_libint_osvrrsxsx_h_

#include <cstdlib>
#include <libint2.h>
#include <util_types.h>
#include <cgshell_ordering.h>

#ifdef __GNUC__
#pragma implementation
#endif

namespace libint2 {

  template <int part, FunctionPosition where, int Lb, int Ld, bool vectorize> struct OSVRR_sx_sx {
    static void compute(const Libint_t* inteval,
        LIBINT2_REALTYPE* target,
        const LIBINT2_REALTYPE* src0,
        const LIBINT2_REALTYPE* src1,
        const LIBINT2_REALTYPE* src2,
        const LIBINT2_REALTYPE* src3,
        const LIBINT2_REALTYPE* src4);
  };

  /** builds (0b|0d)^(m)
      src0 = (0b|0d-1)^(m)
      src1 = (0b|0d-1)^(m+1)
      src2 = (0b|0d-2)^(m)
      src3 = (0b|0d-2)^(m+1)
      src4 = (0b-1|0d-1)^(m+1)
   **/
  template <int Lb, int Ld, bool vectorize> struct OSVRR_sx_sx<0,InKet,Lb,Ld,vectorize> {

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
        const double *QD, *WQ;
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
          inteval->nflops += 8 * Nb * veclen;
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
          inteval->nflops += 3 * Nb * veclen;
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
            inteval->nflops += 3 * veclen;
#endif
            src4_ptr += Ndm1v;

          END_FOR_CART // end of loop over b
        }

        ++id;

      END_FOR_CART // end of loop over d

    }

  };

};

#endif // header guard

