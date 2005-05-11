
#include "libint2_params.h"

#ifndef _libint2_header_
#define _libint2_header_

#define REALTYPE double
#define LIBINT2_MAX_NTARGETS 4

#define VECLEN LIBINT2_MAX_VECLEN

typedef struct {
  REALTYPE __ss_1_over_r_12_ss___up_0[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_1[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_2[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_3[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_4[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_5[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_6[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_7[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_8[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_9[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_10[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_11[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_12[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_13[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_14[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_15[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_16[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_17[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_18[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_19[VECLEN];
  REALTYPE __ss_1_over_r_12_ss___up_20[VECLEN];

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
  
  REALTYPE oo2z[VECLEN];
  REALTYPE oo2e[VECLEN];
  REALTYPE oo2ze[VECLEN];
  REALTYPE roz[VECLEN];
  REALTYPE roe[VECLEN];

  REALTYPE stack[LIBINT2_MAX_STACK_SIZE];
  REALTYPE* targets[LIBINT2_MAX_NTARGETS];

  unsigned int veclength;

} Libint_t;

#endif

#include "libint2_iface.h"

