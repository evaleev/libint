#include <libint2.h>
void compute__ss_1_over_r_12_ds___up_0(Libint_t* libint) {

REALTYPE fp3 = libint->roe * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp2 = libint->__ss_1_over_r_12_ss___up_0 - fp3;
REALTYPE fp1 = libint->oo2e * fp2;
REALTYPE fp41 = libint->WQ_z * libint->__ss_1_over_r_12_ss___up_2;
REALTYPE fp40 = libint->QC_z * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp39 = fp40 + fp41;
REALTYPE fp35 = libint->WQ_z * fp39;
REALTYPE fp38 = libint->WQ_z * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp37 = libint->QC_z * libint->__ss_1_over_r_12_ss___up_0;
REALTYPE fp36 = fp37 + fp38;
REALTYPE fp34 = libint->QC_z * fp36;
REALTYPE fp33 = fp34 + fp35;
REALTYPE fp32 = fp33 + fp1;
libint->stack[5] = fp32;
REALTYPE fp28 = libint->WQ_y * libint->__ss_1_over_r_12_ss___up_2;
REALTYPE fp27 = libint->QC_y * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp26 = fp27 + fp28;
REALTYPE fp31 = libint->WQ_z * fp26;
REALTYPE fp25 = libint->WQ_y * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp24 = libint->QC_y * libint->__ss_1_over_r_12_ss___up_0;
REALTYPE fp23 = fp24 + fp25;
REALTYPE fp30 = libint->QC_z * fp23;
REALTYPE fp29 = fp30 + fp31;
libint->stack[4] = fp29;
REALTYPE fp22 = libint->WQ_y * fp26;
REALTYPE fp21 = libint->QC_y * fp23;
REALTYPE fp20 = fp21 + fp22;
REALTYPE fp19 = fp20 + fp1;
libint->stack[3] = fp19;
REALTYPE fp12 = libint->WQ_x * libint->__ss_1_over_r_12_ss___up_2;
REALTYPE fp11 = libint->QC_x * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp10 = fp11 + fp12;
REALTYPE fp18 = libint->WQ_z * fp10;
REALTYPE fp9 = libint->WQ_x * libint->__ss_1_over_r_12_ss___up_1;
REALTYPE fp8 = libint->QC_x * libint->__ss_1_over_r_12_ss___up_0;
REALTYPE fp7 = fp8 + fp9;
REALTYPE fp17 = libint->QC_z * fp7;
REALTYPE fp16 = fp17 + fp18;
libint->stack[2] = fp16;
REALTYPE fp15 = libint->WQ_y * fp10;
REALTYPE fp14 = libint->QC_y * fp7;
REALTYPE fp13 = fp14 + fp15;
libint->stack[1] = fp13;
REALTYPE fp6 = libint->WQ_x * fp10;
REALTYPE fp5 = libint->QC_x * fp7;
REALTYPE fp4 = fp5 + fp6;
REALTYPE fp0 = fp4 + fp1;
libint->stack[0] = fp0;
/// Number of flops = 42
libint->targets[0] = &(libint->stack[0]);
}

