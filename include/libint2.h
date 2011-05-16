
#ifndef _libint2_header_
#define _libint2_header_

#define LIBINT_T_SS_EREP_SS(mValue) _aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_##mValue
#define LIBINT_T_SS_Km1G12_SS(mValue) _aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_##mValue
#define LIBINT_T_SS_K0G12_SS_0 _aB_s___0__s___1___r12_0_g12_s___0__s___1___Ab__up_0
#define LIBINT_T_SS_K2G12_SS_0 _aB_s___0__s___1___r12_2_g12_s___0__s___1___Ab__up_0
#define LIBINT_T_SS_K4G12_SS_0 _aB_s___0__s___1___r12_4_g12_s___0__s___1___Ab__up_0

#include <libint2_intrinsic_types.h>
#include <libint2_params.h>
#include <libint2_types.h>

#if 0
#define VECLEN LIBINT2_MAX_VECLEN

/** Libint_t is the integrals evaluator object. Libint's evaluator functions take
    pointer to Libint_t as their first argument. The evaluator functions are not reentrant,
    thus each thread should have its own evaluator object. */

typedef struct {
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(0)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(1)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(2)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(3)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(4)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(5)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(6)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(7)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(8)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(9)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(10)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(11)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(12)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(13)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(14)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(15)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(16)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(17)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(18)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(19)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_EREP_SS(20)[VECLEN];

#if SUPPORT_G12
  /**
     Prefactors for recurrence relations from Weber and Daul, Comp. Phys. Comm. 158, 1 (2004).
  */

  LIBINT2_REALTYPE LIBINT_T_SS_K0G12_SS_0[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_K2G12_SS_0[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(0)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(1)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(2)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(3)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(4)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(5)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(6)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(7)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(8)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(9)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(10)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(11)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(12)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(13)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(14)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(15)[VECLEN];
  LIBINT2_REALTYPE LIBINT_T_SS_Km1G12_SS(16)[VECLEN];

  /** LRL1991, Eq. 30, prefactor in front of (a0|c0) */
  LIBINT2_REALTYPE TwoPRepITR_pfac0_0_x[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac0_0_y[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac0_0_z[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac0_1_x[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac0_1_y[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac0_1_z[VECLEN];
  /** LRL1991, Eq. 30, prefactor in front of (a0|c+10) */
  LIBINT2_REALTYPE TwoPRepITR_pfac1_0[VECLEN];
  LIBINT2_REALTYPE TwoPRepITR_pfac1_1[VECLEN];
  /** WD2004, Eq. 30, prefactor in front of (a0|k|c0) */
  LIBINT2_REALTYPE R12kG12_pfac0_0_x[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac0_0_y[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac0_0_z[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac0_1_x[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac0_1_y[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac0_1_z[VECLEN];
  /** WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0) */
  LIBINT2_REALTYPE R12kG12_pfac1_0[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac1_1[VECLEN];
  /** WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0) */
  LIBINT2_REALTYPE R12kG12_pfac2[VECLEN];
  /** WD2004, Eq. 30, prefactor in front of curly brakets (excludes k) */
  LIBINT2_REALTYPE R12kG12_pfac3_0[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac3_1[VECLEN];
  /** WD2004, Eq. 30, prefactor in front of (a0|k-2|c0) */
  LIBINT2_REALTYPE R12kG12_pfac4_0_x[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac4_0_y[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac4_0_z[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac4_1_x[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac4_1_y[VECLEN];
  LIBINT2_REALTYPE R12kG12_pfac4_1_z[VECLEN];

  /** Exponents */
  LIBINT2_REALTYPE zeta_A[VECLEN];
  LIBINT2_REALTYPE zeta_B[VECLEN];
  LIBINT2_REALTYPE zeta_C[VECLEN];
  LIBINT2_REALTYPE zeta_D[VECLEN];
  /** Squared exponents */
  LIBINT2_REALTYPE zeta_A_2[VECLEN];
  LIBINT2_REALTYPE zeta_B_2[VECLEN];
  LIBINT2_REALTYPE zeta_C_2[VECLEN];
  LIBINT2_REALTYPE zeta_D_2[VECLEN];
#endif /* SUPPORT_G12 */

  /**
     Appear in OS RR for ERIs
   */
  /** One over 2.0*zeta */
  LIBINT2_REALTYPE oo2z[VECLEN];
  /** One over 2.0*eta */
  LIBINT2_REALTYPE oo2e[VECLEN];
  /** One over 2.0*(zeta+eta) */
  LIBINT2_REALTYPE oo2ze[VECLEN];
  /** rho over zeta */
  LIBINT2_REALTYPE roz[VECLEN];
  /** rho over eta */
  LIBINT2_REALTYPE roe[VECLEN];

  /** Appear in standard OS RR for ERI and almost all other recurrence relations */
  LIBINT2_REALTYPE WP_x[VECLEN], WP_y[VECLEN], WP_z[VECLEN];
  LIBINT2_REALTYPE WQ_x[VECLEN], WQ_y[VECLEN], WQ_z[VECLEN];
  LIBINT2_REALTYPE PA_x[VECLEN], PA_y[VECLEN], PA_z[VECLEN];
  LIBINT2_REALTYPE QC_x[VECLEN], QC_y[VECLEN], QC_z[VECLEN];
  LIBINT2_REALTYPE AB_x[VECLEN], AB_y[VECLEN], AB_z[VECLEN];
  LIBINT2_REALTYPE CD_x[VECLEN], CD_y[VECLEN], CD_z[VECLEN];
  
  
  /** Scratch buffer to hold intermediates */
  LIBINT2_REALTYPE *stack;
  /** Buffer to hold vector intermediates. Only used by set-level RR
      code if it is vectorized linewise
  */
  LIBINT2_REALTYPE *vstack;
  /** On completion, this contains pointers to computed targets.
      It's modified when library is called, hence should be mutable. */
  mutable LIBINT2_REALTYPE* targets[LIBINT2_MAX_NTARGETS];
  
  /** Actual vector length. Not to exceed VECLEN! If VECLEN is 1 then
      veclength is not used */
  unsigned int veclength;

#if LIBINT2_FLOP_COUNT
  /** FLOP counter. Libint must be configured with --enable-flop-counter
      to allow FLOP counting. It is user's reponsibility to set zero nflops before
      computing integrals. */
  LIBINT2_UINT_LEAST64 nflops;
#endif /* LIBINT2_FLOP_COUNT */

#if LIBINT2_ACCUM_INTS
  /**
     If libint was configured with --enable-accum-ints then the target integrals are accumulated.
     To zero out.the targets automatically before the computation, set this to nonzero.
  */
  int zero_out_targets;
#endif /* LIBINT2_ACCUM_INTS */
  
} Libint_t;
#endif

/**
   these macros define bzero, copy, and inc operations:
*/
/** X[i] = 0 */
#define _libint2_static_api_bzero_short_(X,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] = 0; }
/** X[i] = Y[i] */
#define _libint2_static_api_copy_short_(X,Y,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] = (Y)[i]; }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_short_(X,Y,nelem,a) for(int i=0; i < (nelem); ++i) { (X)[i] = (a) * (Y)[i]; }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_vec_short_(X,Y,nelem,a,vl) for(int i=0, iv=0; i < (nelem)/(vl); ++i) { for(int v=0; v < (vl); ++v, ++iv) { (X)[iv] = (a[v]) * (Y)[iv]; } }
/** X[i] += a*Y[i] */
#define _libint2_static_api_inc_short_(X,Y,nelem,a) for(int i=0; i < (nelem); ++i) { (X)[i] += (a) * (Y)[i]; }
/** X[i] += Y[i] */
#define _libint2_static_api_inc1_short_(X,Y,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] += (Y)[i]; }

#endif

#include "libint2_iface.h"

