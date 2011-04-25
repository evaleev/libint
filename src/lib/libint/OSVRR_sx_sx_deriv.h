
#ifndef _libint2_src_lib_libint_osvrrsxsxderiv_h_
#define _libint2_src_lib_libint_osvrrsxsxderiv_h_

#include <cstdlib>
#include <libint2.h>
#include <util_types.h>
#include <cgshell_ordering.h>

namespace libint2 {

  template <int part, int Lb, int Ld,
            int Da_x,
            int Da_y,
            int Da_z,
            int Db_x,
            int Db_y,
            int Db_z,
            int Dc_x,
            int Dc_y,
            int Dc_z,
            int Dd_x,
            int Dd_y,
            int Dd_z,
            bool vectorize> struct OSVRR_sx_sx_deriv {
    static void compute(const Libint_t* inteval,
        LIBINT2_REALTYPE* target,
        const LIBINT2_REALTYPE* src0,
        const LIBINT2_REALTYPE* src1,
        const LIBINT2_REALTYPE* src2,
        const LIBINT2_REALTYPE* src3,
        const LIBINT2_REALTYPE* src4,
        const LIBINT2_REALTYPE* src5,
        const LIBINT2_REALTYPE* src6,
        const LIBINT2_REALTYPE* src7,
        const LIBINT2_REALTYPE* src8,
        const LIBINT2_REALTYPE* src9,
        const LIBINT2_REALTYPE* src10,
        const LIBINT2_REALTYPE* src11,
        const LIBINT2_REALTYPE* src12,
        const LIBINT2_REALTYPE* src13,
        const LIBINT2_REALTYPE* src14,
        const LIBINT2_REALTYPE* src15,
        const LIBINT2_REALTYPE* src16,
        const LIBINT2_REALTYPE* src17,
        const LIBINT2_REALTYPE* src18,
        const LIBINT2_REALTYPE* src19,
        const LIBINT2_REALTYPE* src20,
        const LIBINT2_REALTYPE* src21,
        const LIBINT2_REALTYPE* src22
        );
  };

  /** builds (a 0|c0)^(m)
      src0 = (a-10|c0)^(m)
      src1 = (a-10|c0)^(m+1)
      src2 = (a-20|c0)^(m)
      src3 = (a-20|c0)^(m+1)
      src4 = (a-10|c-10)^(m+1)

      src5 = Da_x-1 (a-10|c0)^(m)
      src6 = Da_x-1 (a-10|c0)^(m+1)
      src7 = Da_y-1 (a-10|c0)^(m)
      src8 = Da_y-1 (a-10|c0)^(m+1)
      src9 = Da_z-1 (a-10|c0)^(m)
      src10= Da_z-1 (a-10|c0)^(m+1)

      src11= Db_x-1 (a-10|c0)^(m)
      src12= Db_x-1 (a-10|c0)^(m+1)
      src13= Db_y-1 (a-10|c0)^(m)
      src14= Db_y-1 (a-10|c0)^(m+1)
      src15= Db_z-1 (a-10|c0)^(m)
      src16= Db_z-1 (a-10|c0)^(m+1)

      src17= Dc_x-1 (a-10|c0)^(m+1)
      src18= Dc_y-1 (a-10|c0)^(m+1)
      src19= Dc_z-1 (a-10|c0)^(m+1)

      src20= Dd_x-1 (a-10|c0)^(m+1)
      src21= Dd_y-1 (a-10|c0)^(m+1)
      src22= Dd_z-1 (a-10|c0)^(m+1)
   **/
  template <int Lb, int Ld,
  int Da_x,
  int Da_y,
  int Da_z,
  int Db_x,
  int Db_y,
  int Db_z,
  int Dc_x,
  int Dc_y,
  int Dc_z,
  int Dd_x,
  int Dd_y,
  int Dd_z,
  bool vectorize> struct OSVRR_sx_sx_deriv<0,Lb,Ld,
                                     Da_x,Da_y,Da_z,
                                     Db_x,Db_y,Db_z,
                                     Dc_x,Dc_y,Dc_z,
                                     Dd_x,Dd_y,Dd_z,
                                     vectorize> {

    static void compute(const Libint_t* inteval,
        LIBINT2_REALTYPE* target,
        const LIBINT2_REALTYPE* src0,
        const LIBINT2_REALTYPE* src1,
        const LIBINT2_REALTYPE* src2,
        const LIBINT2_REALTYPE* src3,
        const LIBINT2_REALTYPE* src4,
        const LIBINT2_REALTYPE* src5,
        const LIBINT2_REALTYPE* src6,
        const LIBINT2_REALTYPE* src7,
        const LIBINT2_REALTYPE* src8,
        const LIBINT2_REALTYPE* src9,
        const LIBINT2_REALTYPE* src10,
        const LIBINT2_REALTYPE* src11,
        const LIBINT2_REALTYPE* src12,
        const LIBINT2_REALTYPE* src13,
        const LIBINT2_REALTYPE* src14,
        const LIBINT2_REALTYPE* src15,
        const LIBINT2_REALTYPE* src16,
        const LIBINT2_REALTYPE* src17,
        const LIBINT2_REALTYPE* src18,
        const LIBINT2_REALTYPE* src19,
        const LIBINT2_REALTYPE* src20,
        const LIBINT2_REALTYPE* src21,
        const LIBINT2_REALTYPE* src22
        ) {

      // works for (sd|sp) and higher
      if (Lb < 2 || Ld < 1)
        abort();

      const unsigned int veclen = vectorize ? inteval->veclen : 1;

      const unsigned int Nd = INT_NCART(Ld);
      const unsigned int NdV = Nd * veclen;

      int bx, by, bz;
      FOR_CART(bx, by, bz, Lb)

        int b[3]; b[0] = bx;  b[1] = by;  b[2] = bz;

        enum XYZ {x=0, y=1, z=2};
        // Build along x, if possible
        XYZ xyz = z;
        if (by != 0) xyz = y;
        if (bx != 0) xyz = x;
        --b[xyz];

        // redirect
        const double *PB, *WP;
        switch(xyz) {
          case x:
            PB = inteval->PB_x;
            WP = inteval->WP_x;
            break;
          case y:
            PB = inteval->PB_y;
            WP = inteval->WP_y;
            break;
          case z:
            PB = inteval->PB_z;
            WP = inteval->WP_z;
            break;
        }

        const unsigned int ibm1 = INT_CARTINDEX(Lb-1,b[0],b[1]);
        const unsigned int bm10d0_offset = ibm1 * NdV;
        const LIBINT2_REALTYPE* src0_ptr = src0 + bm10d0_offset;
        const LIBINT2_REALTYPE* src1_ptr = src1 + bm10d0_offset;

        // if b-2_xyz exists, include (0 b-2_xyz | 0 d)
        if (b[xyz] > 0) {
          --b[xyz];
          const unsigned int ibm2 = INT_CARTINDEX(Lb-2,b[0],b[1]);
          const unsigned int bm20d0_offset = ibm2 * NdV;
          ++b[xyz];
          const LIBINT2_REALTYPE* src2_ptr = src2 + bm20d0_offset;
          const LIBINT2_REALTYPE* src3_ptr = src3 + bm20d0_offset;
          const LIBINT2_REALTYPE bxyz = (LIBINT2_REALTYPE)b[xyz];

          unsigned int dv = 0;
          for(unsigned int d = 0; d < Nd; ++d) {
            for(unsigned int v=0; v<veclen; ++v, ++dv) {
              target[dv] = PB[v] * src0_ptr[dv] + WP[v] * src1_ptr[dv] + bxyz * inteval->oo2z[v] * (src2_ptr[dv] - inteval->roz[v] * src3_ptr[dv]);
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 8 * NdV;
#endif

        }
        else {
          unsigned int dv = 0;
          for(unsigned int d = 0; d < Nd; ++d) {
            for(unsigned int v=0; v<veclen; ++v, ++dv) {
              target[dv] = PB[v] * src0_ptr[dv] + WP[v] * src1_ptr[dv];
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }

        {
          const unsigned int Ndm1 = INT_NCART(Ld-1);
          const unsigned int Ndm1V = Ndm1 * veclen;
          const unsigned int bm10dm10_offset = ibm1 * Ndm1V;
          const LIBINT2_REALTYPE* src4_ptr = src4 + bm10dm10_offset;

          // loop over d-1 shell and include (0 b-1_xyz | 0 d-1_xyz) to (0 b | 0 d)
          int dx, dy, dz;
          FOR_CART(dx, dy, dz, Ld-1)

            int d[3]; d[0] = dx;  d[1] = dy;  d[2] = dz;
            ++d[xyz];

            const unsigned int dc = INT_CARTINDEX(Ld,d[0],d[1]);
            const unsigned int dc_offset = dc * veclen;
            LIBINT2_REALTYPE* tptr = target + dc_offset;
            const LIBINT2_REALTYPE dxyz = (LIBINT2_REALTYPE)d[xyz];
            for(unsigned int v=0; v<veclen; ++v) {
              tptr[v] += dxyz * inteval->oo2ze[v] * src4_ptr[v];
            }
#if LIBINT2_FLOP_COUNT
            inteval->nflops += 3 * veclen;
#endif
            src4_ptr += veclen;

          END_FOR_CART
        }

#define dcontrA(target,srcA,srcB,id,ecoef1,ecoef2) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + bm10d0_offset; \
  const LIBINT2_REALTYPE* srcB_ptr = srcB + bm10d0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  const LIBINT2_REALTYPE* c2 = inteval->ecoef2; \
  unsigned int dv = 0; \
  for(unsigned int d = 0; d < Nd; ++d) { \
    for(unsigned int v=0; v<veclen; ++v, ++dv) { \
      target[dv] +=  di * (c1[v] * srcA_ptr[dv] - c2[v] * srcB_ptr[dv]); \
    } \
  } \
}

#define dcontrB(target,srcA,srcB,id,ecoef1,ecoef2) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + bm10d0_offset; \
  const LIBINT2_REALTYPE* srcB_ptr = srcB + bm10d0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  const LIBINT2_REALTYPE* c2 = inteval->ecoef2; \
  unsigned int dv = 0; \
  for(unsigned int d = 0; d < Nd; ++d) { \
    for(unsigned int v=0; v<veclen; ++v, ++dv) { \
      target[dv] -=  di * (c1[v] * srcA_ptr[dv] + c2[v] * srcB_ptr[dv]); \
    } \
  } \
}

#define dcontrCD(target,srcA,id,ecoef1) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + bm10d0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  unsigned int dv = 0; \
  for(unsigned int d = 0; d < Nd; ++d) { \
    for(unsigned int v=0; v<veclen; ++v, ++dv) { \
      target[dv] +=  di * c1[v] * srcA_ptr[dv]; \
    } \
  } \
}

        // if Da_x-1 exists
#if LIBINT2_DEFINED(any,rho12_over_alpha2) && LIBINT2_DEFINED(any,alpha1_rho_over_zeta2)
        if (Da_x > 0 && xyz == x){
          dcontrA(target,src5,src6,Da_x,rho12_over_alpha2,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
        if (Da_y > 0 && xyz == y){
          dcontrA(target,src11,src12,Da_y,rho12_over_alpha2,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
        if (Da_z > 0 && xyz == z){
          dcontrA(target,src17,src18,Da_z,rho12_over_alpha2,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
#endif
        // if Db_x-1 exists
#if LIBINT2_DEFINED(any,rho12_over_alpha2) && LIBINT2_DEFINED(any,alpha2_rho_over_zeta2)
        if (Db_x > 0 && xyz == x){
          dcontrB(target,src7,src8,Db_x,rho12_over_alpha2,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
        if (Db_y > 0 && xyz == y){
          dcontrB(target,src13,src14,Db_y,rho12_over_alpha2,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
        if (Db_z > 0 && xyz == z){
          dcontrB(target,src19,src20,Db_z,rho12_over_alpha2,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NdV;
#endif
        }
#endif
        // if Dc_x-1 exists
#if LIBINT2_DEFINED(any,alpha3_over_zetapluseta)
        if (Dc_x > 0 && xyz == x){
          dcontrCD(target,src9,Dc_x,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
        if (Dc_y > 0 && xyz == y){
          dcontrCD(target,src15,Dc_y,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
        if (Dc_z > 0 && xyz == z){
          dcontrCD(target,src21,Dc_z,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
#endif
        // if Dd_x-1 exists
#if LIBINT2_DEFINED(any,alpha4_over_zetapluseta)
        if (Dd_x > 0 && xyz == x){
          dcontrCD(target,src10,Dd_x,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
        if (Dd_y > 0 && xyz == y){
          dcontrCD(target,src16,Dd_y,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
        if (Dd_z > 0 && xyz == z){
          dcontrCD(target,src22,Dd_z,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NdV;
#endif
        }
#endif

        target += NdV;

      END_FOR_CART // end of loop over a-1

      /** Number of flops = ??? */
      //inteval->nflops = inteval->nflops + 222 * 1 * 1 * veclen;

    }

  };

};

#endif // header guard

