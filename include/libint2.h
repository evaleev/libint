
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

  REALTYPE WP_x[VECLEN], WP_y[VECLEN], WP_z[VECLEN];
  REALTYPE WQ_x[VECLEN], WQ_y[VECLEN], WQ_z[VECLEN];
  REALTYPE PA_x[VECLEN], PA_y[VECLEN], PA_z[VECLEN];
  REALTYPE QC_x[VECLEN], QC_y[VECLEN], QC_z[VECLEN];
  REALTYPE AB_x[VECLEN], AB_y[VECLEN], AB_z[VECLEN];
  REALTYPE CD_x[VECLEN], CD_y[VECLEN], CD_z[VECLEN];

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

