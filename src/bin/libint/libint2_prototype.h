
#ifndef _libint2_header_
#define _libint2_header_

#define REALTYPE double
#define LIBINT_MAX_STACK 1000000
#define LIBINT_MAX_NTARGETS 4

typedef struct {
  REALTYPE __ss_1_over_r_12_ss___up_0;
  REALTYPE __ss_1_over_r_12_ss___up_1;
  REALTYPE __ss_1_over_r_12_ss___up_2;
  REALTYPE __ss_1_over_r_12_ss___up_3;
  REALTYPE __ss_1_over_r_12_ss___up_4;
  REALTYPE __ss_1_over_r_12_ss___up_5;
  REALTYPE __ss_1_over_r_12_ss___up_6;

  REALTYPE WP_x, WP_y, WP_z;
  REALTYPE WQ_x, WQ_y, WQ_z;
  REALTYPE PA_x, PA_y, PA_z;
  REALTYPE QC_x, QC_y, QC_z;
  REALTYPE AB_x, AB_y, AB_z;
  REALTYPE CD_x, CD_y, CD_z;

  REALTYPE oo2z;
  REALTYPE oo2e;
  REALTYPE oo2ze;
  REALTYPE roz;
  REALTYPE roe;

  REALTYPE* stack;
  REALTYPE* targets[LIBINT_MAX_NTARGETS];

} Libint_t;

#endif

