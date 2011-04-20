
#ifndef _libint2_src_lib_libint_osvrrxsxsderiv_h_
#define _libint2_src_lib_libint_osvrrxsxsderiv_h_

#include <cstdlib>
#include <libint2.h>
#include <util_types.h>
#include <cgshell_ordering.h>

namespace libint2 {

  template <int part, int La, int Lc,
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
            bool vectorize> struct OSVRR_xs_xs_deriv {
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
  template <int La, int Lc,
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
  bool vectorize> struct OSVRR_xs_xs_deriv<0,La,Lc,
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

      // works for (ds|ps) and higher
      if (La < 2 || Lc < 1)
        abort();

      const unsigned int veclen = vectorize ? inteval->veclen : 1;

      const unsigned int Nc = INT_NCART(Lc);
      const unsigned int NcV = Nc * veclen;

      int ax, ay, az;
      FOR_CART(ax, ay, az, La)

        int a[3]; a[0] = ax;  a[1] = ay;  a[2] = az;

        enum XYZ {x=0, y=1, z=2};
        // Build along x, if possible
        XYZ xyz = z;
        if (ay != 0) xyz = y;
        if (ax != 0) xyz = x;
        --a[xyz];

        // redirect
        const double *PA, *WP;
        switch(xyz) {
          case x:
            PA = inteval->PA_x;
            WP = inteval->WP_x;
            break;
          case y:
            PA = inteval->PA_y;
            WP = inteval->WP_y;
            break;
          case z:
            PA = inteval->PA_z;
            WP = inteval->WP_z;
            break;
        }

        const unsigned int iam1 = INT_CARTINDEX(La-1,a[0],a[1]);
        const unsigned int am10c0_offset = iam1 * NcV;
        const LIBINT2_REALTYPE* src0_ptr = src0 + am10c0_offset;
        const LIBINT2_REALTYPE* src1_ptr = src1 + am10c0_offset;

        // if a-2_xyz exists, include (a-2_xyz 0 | c 0)
        if (a[xyz] > 0) {
          --a[xyz];
          const unsigned int iam2 = INT_CARTINDEX(La-2,a[0],a[1]);
          const unsigned int am20c0_offset = iam2 * NcV;
          ++a[xyz];
          const LIBINT2_REALTYPE* src2_ptr = src2 + am20c0_offset;
          const LIBINT2_REALTYPE* src3_ptr = src3 + am20c0_offset;
          const LIBINT2_REALTYPE axyz = (LIBINT2_REALTYPE)a[xyz];

          unsigned int cv = 0;
          for(unsigned int c = 0; c < Nc; ++c) {
            for(unsigned int v=0; v<veclen; ++v, ++cv) {
              target[cv] = PA[v] * src0_ptr[cv] + WP[v] * src1_ptr[cv] + axyz * inteval->oo2z[v] * (src2_ptr[cv] - inteval->roz[v] * src3_ptr[cv]);
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 8 * NcV;
#endif

        }
        else {
          unsigned int cv = 0;
          for(unsigned int c = 0; c < Nc; ++c) {
            for(unsigned int v=0; v<veclen; ++v, ++cv) {
              target[cv] = PA[v] * src0_ptr[cv] + WP[v] * src1_ptr[cv];
            }
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }

        {
          const unsigned int Ncm1 = INT_NCART(Lc-1);
          const unsigned int Ncm1V = Ncm1 * veclen;
          const unsigned int am10cm10_offset = iam1 * Ncm1V;
          const LIBINT2_REALTYPE* src4_ptr = src4 + am10cm10_offset;

          // loop over c-1 shell and include (a-1_xyz 0 | c-1_xyz 0) to (a 0 | c 0)
          int cx, cy, cz;
          FOR_CART(cx, cy, cz, Lc-1)

            int c[3]; c[0] = cx;  c[1] = cy;  c[2] = cz;
            ++c[xyz];

            const unsigned int cc = INT_CARTINDEX(Lc,c[0],c[1]);
            const unsigned int cc_offset = cc * veclen;
            LIBINT2_REALTYPE* tptr = target + cc_offset;
            const LIBINT2_REALTYPE cxyz = (LIBINT2_REALTYPE)c[xyz];
            for(unsigned int v=0; v<veclen; ++v) {
              tptr[v] += cxyz * inteval->oo2ze[v] * src4_ptr[v];
            }
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * veclen;
#endif
            src4_ptr += veclen;

          END_FOR_CART // end of loop over c-1
        }


#define dcontrA(target,srcA,srcB,id,ecoef1,ecoef2) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + am10c0_offset; \
  const LIBINT2_REALTYPE* srcB_ptr = srcB + am10c0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  const LIBINT2_REALTYPE* c2 = inteval->ecoef2; \
  unsigned int cv = 0; \
  for(unsigned int c = 0; c < Nc; ++c) { \
    for(unsigned int v=0; v<veclen; ++v, ++cv) { \
      target[cv] -=  di * (c1[v] * srcA_ptr[cv] + c2[v] * srcB_ptr[cv]); \
    } \
  } \
}

#define dcontrB(target,srcA,srcB,id,ecoef1,ecoef2) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + am10c0_offset; \
  const LIBINT2_REALTYPE* srcB_ptr = srcB + am10c0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  const LIBINT2_REALTYPE* c2 = inteval->ecoef2; \
  unsigned int cv = 0; \
  for(unsigned int c = 0; c < Nc; ++c) { \
    for(unsigned int v=0; v<veclen; ++v, ++cv) { \
      target[cv] +=  di * (c1[v] * srcA_ptr[cv] - c2[v] * srcB_ptr[cv]); \
    } \
  } \
}

#define dcontrCD(target,srcA,id,ecoef1) { \
  const LIBINT2_REALTYPE* srcA_ptr = srcA + am10c0_offset; \
  const LIBINT2_REALTYPE di = (LIBINT2_REALTYPE)id; \
  const LIBINT2_REALTYPE* c1 = inteval->ecoef1; \
  unsigned int cv = 0; \
  for(unsigned int c = 0; c < Nc; ++c) { \
    for(unsigned int v=0; v<veclen; ++v, ++cv) { \
      target[cv] +=  di * c1[v] * srcA_ptr[cv]; \
    } \
  } \
}

        // if Da_x-1 exists
#if LIBINT2_DEFINED(any,rho12_over_alpha1) && LIBINT2_DEFINED(any,alpha1_rho_over_zeta2)
        if (Da_x > 0 && xyz == x){
          dcontrA(target,src5,src6,Da_x,rho12_over_alpha1,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
        if (Da_y > 0 && xyz == y){
          dcontrA(target,src7,src8,Da_y,rho12_over_alpha1,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
        if (Da_z > 0 && xyz == z){
          dcontrA(target,src9,src10,Da_z,rho12_over_alpha1,alpha1_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
#endif
        // if Db_x-1 exists
#if LIBINT2_DEFINED(any,rho12_over_alpha1) && LIBINT2_DEFINED(any,alpha2_rho_over_zeta2)
        if (Db_x > 0 && xyz == x){
          dcontrB(target,src11,src12,Db_x,rho12_over_alpha1,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
        if (Db_y > 0 && xyz == y){
          dcontrB(target,src13,src14,Db_y,rho12_over_alpha1,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
        if (Db_z > 0 && xyz == z){
          dcontrB(target,src15,src16,Db_z,rho12_over_alpha1,alpha2_rho_over_zeta2);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 5 * NcV;
#endif
        }
#endif
        // if Dc_x-1 exists
#if LIBINT2_DEFINED(any,alpha3_over_zetapluseta)
        if (Dc_x > 0 && xyz == x){
          dcontrCD(target,src17,Dc_x,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
        if (Dc_y > 0 && xyz == y){
          dcontrCD(target,src18,Dc_y,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
        if (Dc_z > 0 && xyz == z){
          dcontrCD(target,src19,Dc_z,alpha3_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
#endif
        // if Dd_x-1 exists
#if LIBINT2_DEFINED(any,alpha4_over_zetapluseta)
        if (Dd_x > 0 && xyz == x){
          dcontrCD(target,src20,Dd_x,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
        if (Dd_y > 0 && xyz == y){
          dcontrCD(target,src21,Dd_y,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
        if (Dd_z > 0 && xyz == z){
          dcontrCD(target,src22,Dd_z,alpha4_over_zetapluseta);
#if LIBINT2_FLOP_COUNT
          inteval->nflops += 3 * NcV;
#endif
        }
#endif

        target += NcV;

      END_FOR_CART // end of loop over a-1

      /** Number of flops = ??? */
      //inteval->nflops = inteval->nflops + 222 * 1 * 1 * veclen;

    }

  };

};

#endif // header guard

