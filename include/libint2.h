
#include <libint2_types.h>
#include <libint2_params.h>

#ifndef _libint2_header_
#define _libint2_header_

#define REALTYPE double
#define LIBINT2_MAX_NTARGETS 4

#define VECLEN LIBINT2_MAX_VECLEN

#define  LIBINT_T_SS_EREP_SS(mValue) _aB_s__0__s__1___TwoERep_s__0__s__1___Ab__up_##mValue

/** Libint_t is the integrals evaluator object. Libint's evaluator functions take
    pointer to Libint_t as their first argument. The evaluator functions are not reentrant,
    thus each thread should have its own evaluator object. */

typedef struct {
  REALTYPE LIBINT_T_SS_EREP_SS(0)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(1)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(2)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(3)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(4)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(5)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(6)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(7)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(8)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(9)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(10)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(11)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(12)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(13)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(14)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(15)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(16)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(17)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(18)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(19)[VECLEN];
  REALTYPE LIBINT_T_SS_EREP_SS(20)[VECLEN];

  //
  // Prefactors for recurrence relations from Weber and Daul, Comp. Phys. Comm. 158, 1 (2004).
  //

  REALTYPE __ss_r_12_up_2_times_G12_ss___up_0[VECLEN];
  REALTYPE __ss_r_12_up_0_times_G12_ss___up_0[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_0[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_1[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_2[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_3[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_4[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_5[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_6[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_7[VECLEN];
  REALTYPE __ss_r_12_up__minus_1_times_G12_ss___up_8[VECLEN];
  /// WD2004, Eq. 30, prefactor in front of (a0|k|c0)
  REALTYPE R12kG12_pfac0_0_x[VECLEN];
  REALTYPE R12kG12_pfac0_0_y[VECLEN];
  REALTYPE R12kG12_pfac0_0_z[VECLEN];
  REALTYPE R12kG12_pfac0_1_x[VECLEN];
  REALTYPE R12kG12_pfac0_1_y[VECLEN];
  REALTYPE R12kG12_pfac0_1_z[VECLEN];
  /// WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0)
  REALTYPE R12kG12_pfac1_0[VECLEN];
  REALTYPE R12kG12_pfac1_1[VECLEN];
  /// WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0)
  REALTYPE R12kG12_pfac2[VECLEN];
  /// WD2004, Eq. 30, prefactor in front of curly brakets (excludes k)
  REALTYPE R12kG12_pfac3_0[VECLEN];
  REALTYPE R12kG12_pfac3_1[VECLEN];
  /// WD2004, Eq. 30, prefactor in front of (a0|k-2|c0)
  REALTYPE R12kG12_pfac4_0_x[VECLEN];
  REALTYPE R12kG12_pfac4_0_y[VECLEN];
  REALTYPE R12kG12_pfac4_0_z[VECLEN];
  REALTYPE R12kG12_pfac4_1_x[VECLEN];
  REALTYPE R12kG12_pfac4_1_y[VECLEN];
  REALTYPE R12kG12_pfac4_1_z[VECLEN];
  
  REALTYPE WP_x[VECLEN], WP_y[VECLEN], WP_z[VECLEN];
  REALTYPE WQ_x[VECLEN], WQ_y[VECLEN], WQ_z[VECLEN];
  REALTYPE PA_x[VECLEN], PA_y[VECLEN], PA_z[VECLEN];
  REALTYPE QC_x[VECLEN], QC_y[VECLEN], QC_z[VECLEN];
  REALTYPE AB_x[VECLEN], AB_y[VECLEN], AB_z[VECLEN];
  REALTYPE CD_x[VECLEN], CD_y[VECLEN], CD_z[VECLEN];
  
  /// Exponents
  REALTYPE zeta_A[VECLEN];
  REALTYPE zeta_B[VECLEN];
  REALTYPE zeta_C[VECLEN];
  REALTYPE zeta_D[VECLEN];
  /// Squared exponents
  REALTYPE zeta_A_2[VECLEN];
  REALTYPE zeta_B_2[VECLEN];
  REALTYPE zeta_C_2[VECLEN];
  REALTYPE zeta_D_2[VECLEN];
  
  /// One over 2.0*zeta
  REALTYPE oo2z[VECLEN];
  /// One over 2.0*eta
  REALTYPE oo2e[VECLEN];
  /// One over 2.0*(zeta+eta)
  REALTYPE oo2ze[VECLEN];
  /// rho over zeta
  REALTYPE roz[VECLEN];
  /// rho over eta
  REALTYPE roe[VECLEN];
  
  /// Stack of the intermediates is here
  REALTYPE stack[LIBINT2_MAX_STACK_SIZE * VECLEN];
  /// On completion, this contains pointers to computed targets
  REALTYPE* targets[LIBINT2_MAX_NTARGETS];
  
  /** Actual vector length. Not to exceed VECLEN! If VECLEN is 1 then
      veclength is not used */
  unsigned int veclength;
  /** FLOP counter. Libint must be configured with --enable-flop-counter
      to allow FLOP counting. It is user's reponsibility to set zero nflops before
      computing integrals. */
  LIBINT2_UINT_LEAST64 nflops;
  
} Libint_t;

#endif

#include "libint2_iface.h"

