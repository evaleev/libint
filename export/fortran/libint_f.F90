MODULE libint_f
   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_PTR, C_NULL_PTR, C_INT, C_FUNPTR, C_F_PROCPOINTER, C_SIZE_T

#include <config.h>
#include <libint2_params.h>
#include "fortran_incldefs.h"

   IMPLICIT NONE

#ifdef LIBINT2_MAX_AM
   INTEGER, PARAMETER :: libint2_max_am = LIBINT2_MAX_AM
#endif
#ifdef LIBINT2_MAX_AM_eri
   INTEGER, PARAMETER :: libint2_max_am_eri = LIBINT2_MAX_AM_eri
#endif
#ifdef LIBINT2_MAX_AM_eri1
   INTEGER, PARAMETER :: libint2_max_am_eri1 = LIBINT2_MAX_AM_eri1
#endif
#ifdef LIBINT2_MAX_AM_eri2
   INTEGER, PARAMETER :: libint2_max_am_eri2 = LIBINT2_MAX_AM_eri2
#endif

   INTEGER, PARAMETER :: libint2_max_veclen = LIBINT2_MAX_VECLEN

   TYPE, BIND(C) :: libint_t
#ifdef LIBINT2_DEFINED_PO2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PO2
#endif
#ifdef LIBINT2_DEFINED_R12_2_G12_scale_to_G12T1G12
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12_2_G12_scale_to_G12T1G12
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac2
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac1_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac1_0
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac1_1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac1_1
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac3_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac3_0
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac3_1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac3_1
#endif

#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_1
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_2
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_3
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_3
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_4
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_4
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_5
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_5
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_6
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_6
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_7
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_7
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_8
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_8
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_9
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_9
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_10
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_10
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_11
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_11
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_12
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_12
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_13
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_13
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_14
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_14
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_15
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_15
#endif
#ifdef LIBINT2_DEFINED__aB_s___0___ElecPot_s___0___Ab__up_16
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0___ElecPot_s___0___Ab__up_16
#endif

#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_21
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_21
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_22
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_22
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_23
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_23
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_24
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_24
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_25
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_25
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_26
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_26
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_27
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_27
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_28
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_28
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_29
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_29
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_30
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_30
#endif

#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_0_g12_s___0__s___1___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_0_g12_s___0__s___1___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_2_g12_s___0__s___1___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_2_g12_s___0__s___1___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_4_g12_s___0__s___1___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_4_g12_s___0__s___1___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_0
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_0
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_1
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_2
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_3
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_3
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_4
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_4
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_5
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_5
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_6
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_6
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_7
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_7
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_8
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_8
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_9
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_9
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_10
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_10
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_11
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_11
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_12
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_12
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_13
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_13
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_14
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_14
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_15
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_15
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_16
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_16
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_17
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_17
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_18
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_18
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_19
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_19
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_20
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_20
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_21
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_21
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_22
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_22
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_23
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_23
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_24
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_24
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_25
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_25
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_26
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_26
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_27
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_27
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_28
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_28
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_29
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_29
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_30
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_30
#endif

#ifdef LIBINT2_DEFINED_alpha1_rho_over_zeta2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha1_rho_over_zeta2
#endif
#ifdef LIBINT2_DEFINED_alpha2_rho_over_zeta2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha2_rho_over_zeta2
#endif
#ifdef LIBINT2_DEFINED_alpha3_rho_over_eta2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha3_rho_over_eta2
#endif
#ifdef LIBINT2_DEFINED_alpha4_rho_over_eta2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha4_rho_over_eta2
#endif
#ifdef LIBINT2_DEFINED_rho12_over_alpha1
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: rho12_over_alpha1
#endif
#ifdef LIBINT2_DEFINED_rho12_over_alpha2
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: rho12_over_alpha2
#endif
#ifdef LIBINT2_DEFINED_rho34_over_alpha3
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: rho34_over_alpha3
#endif
#ifdef LIBINT2_DEFINED_rho34_over_alpha4
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: rho34_over_alpha4
#endif
#ifdef LIBINT2_DEFINED_AB_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AB_x
#endif
#ifdef LIBINT2_DEFINED_AB_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AB_y
#endif
#ifdef LIBINT2_DEFINED_AB_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AB_z
#endif
#ifdef LIBINT2_DEFINED_AC_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AC_x
#endif
#ifdef LIBINT2_DEFINED_AC_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AC_y
#endif
#ifdef LIBINT2_DEFINED_AC_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: AC_z
#endif
#ifdef LIBINT2_DEFINED_BA_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BA_x
#endif
#ifdef LIBINT2_DEFINED_BA_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BA_y
#endif
#ifdef LIBINT2_DEFINED_BA_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BA_z
#endif
#ifdef LIBINT2_DEFINED_BD_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BD_x
#endif
#ifdef LIBINT2_DEFINED_BD_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BD_y
#endif
#ifdef LIBINT2_DEFINED_BD_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BD_z
#endif
#ifdef LIBINT2_DEFINED_BO_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BO_x
#endif
#ifdef LIBINT2_DEFINED_BO_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BO_y
#endif
#ifdef LIBINT2_DEFINED_BO_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: BO_z
#endif
#ifdef LIBINT2_DEFINED_CD_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: CD_x
#endif
#ifdef LIBINT2_DEFINED_CD_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: CD_y
#endif
#ifdef LIBINT2_DEFINED_CD_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: CD_z
#endif
#ifdef LIBINT2_DEFINED_DC_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) DC_x
#endif
#ifdef LIBINT2_DEFINED_DC_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) DC_y
#endif
#ifdef LIBINT2_DEFINED_DC_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) DC_z
#endif
#ifdef LIBINT2_DEFINED_PA_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PA_x
#endif
#ifdef LIBINT2_DEFINED_PA_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PA_y
#endif
#ifdef LIBINT2_DEFINED_PA_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PA_z
#endif
#ifdef LIBINT2_DEFINED_PB_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PB_x
#endif
#ifdef LIBINT2_DEFINED_PB_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PB_y
#endif
#ifdef LIBINT2_DEFINED_PB_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PB_z
#endif
#ifdef LIBINT2_DEFINED_PC_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PC_x
#endif
#ifdef LIBINT2_DEFINED_PC_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PC_y
#endif
#ifdef LIBINT2_DEFINED_PC_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PC_z
#endif
#ifdef LIBINT2_DEFINED_PO_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PO_x
#endif
#ifdef LIBINT2_DEFINED_PO_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PO_y
#endif
#ifdef LIBINT2_DEFINED_PO_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: PO_z
#endif
#ifdef LIBINT2_DEFINED_QC_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QC_x
#endif
#ifdef LIBINT2_DEFINED_QC_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QC_y
#endif
#ifdef LIBINT2_DEFINED_QC_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QC_z
#endif
#ifdef LIBINT2_DEFINED_QD_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QD_x
#endif
#ifdef LIBINT2_DEFINED_QD_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QD_y
#endif
#ifdef LIBINT2_DEFINED_QD_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: QD_z
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_0_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_0_x
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_0_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_0_y
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_0_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_0_z
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_1_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_1_x
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_1_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_1_y
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac0_1_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac0_1_z
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_0_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_0_x
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_0_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_0_y
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_0_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_0_z
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_1_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_1_x
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_1_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_1_y
#endif
#ifdef LIBINT2_DEFINED_R12kG12_pfac4_1_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: R12kG12_pfac4_1_z
#endif
#ifdef LIBINT2_DEFINED_WP_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WP_x
#endif
#ifdef LIBINT2_DEFINED_WP_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WP_y
#endif
#ifdef LIBINT2_DEFINED_WP_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WP_z
#endif
#ifdef LIBINT2_DEFINED_WQ_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WQ_x
#endif
#ifdef LIBINT2_DEFINED_WQ_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WQ_y
#endif
#ifdef LIBINT2_DEFINED_WQ_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: WQ_z
#endif
#ifdef LIBINT2_DEFINED__0_Overlap_0_x
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_0_Overlap_0_x
#endif
#ifdef LIBINT2_DEFINED__0_Overlap_0_y
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_0_Overlap_0_y
#endif
#ifdef LIBINT2_DEFINED__0_Overlap_0_z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: f_0_Overlap_0_z
#endif
#ifdef LIBINT2_DEFINED_alpha1_over_zetapluseta
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha1_over_zetapluseta
#endif
#ifdef LIBINT2_DEFINED_alpha2_over_zetapluseta
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha2_over_zetapluseta
#endif
#ifdef LIBINT2_DEFINED_alpha3_over_zetapluseta
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha3_over_zetapluseta
#endif
#ifdef LIBINT2_DEFINED_alpha4_over_zetapluseta
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: alpha4_over_zetapluseta
#endif
#ifdef LIBINT2_DEFINED_gamma
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: gamma
#endif
#ifdef LIBINT2_DEFINED_gamma_bra
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: gamma_bra
#endif
#ifdef LIBINT2_DEFINED_gamma_ket
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: gamma_ket
#endif
#ifdef LIBINT2_DEFINED_oo2e
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: oo2e
#endif
#ifdef LIBINT2_DEFINED_oo2z
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: oo2z
#endif
#ifdef LIBINT2_DEFINED_oo2ze
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: oo2ze
#endif
#ifdef LIBINT2_DEFINED_roe
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: roe
#endif
#ifdef LIBINT2_DEFINED_roz
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: roz
#endif
#ifdef LIBINT2_DEFINED_two_alpha0_bra
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: two_alpha0_bra
#endif
#ifdef LIBINT2_DEFINED_two_alpha0_ket
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: two_alpha0_ket
#endif
#ifdef LIBINT2_DEFINED_two_alpha1_bra
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: two_alpha1_bra
#endif
#ifdef LIBINT2_DEFINED_two_alpha1_ket
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: two_alpha1_ket
#endif
#ifdef LIBINT2_DEFINED_zeta_A
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: zeta_A
#endif
#ifdef LIBINT2_DEFINED_zeta_B
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: zeta_B
#endif
#ifdef LIBINT2_DEFINED_zeta_C
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: zeta_C
#endif
#ifdef LIBINT2_DEFINED_zeta_D
      REAL(KIND=C_DOUBLE), DIMENSION(libint2_max_veclen) :: zeta_D
#endif
      TYPE(C_PTR) :: stack
      TYPE(C_PTR) :: vstack
      TYPE(C_PTR), DIMENSION(LIBINT2_MAX_NTARGETS) :: targets
      INTEGER(KIND=C_INT) :: veclen
#if LIBINT2_FLOP_COUNT
      TYPE(C_PTR) :: nflops
#endif
#if LIBINT2_PROFILE
#if LIBINT2_CPLUSPLUS_STD >= 2011
      TYPE(C_PTR) :: timers
#endif
#endif
#if LIBINT2_ACCUM_INTS
      INTEGER(KIND=C_INT) :: zero_out_targets
#endif
#if LIBINT2_CONTRACTED_INTS
      INTEGER(KIND=C_INT) :: contrdepth
#endif
   END TYPE

#ifdef INCLUDE_ERI
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri), &
      BIND(C) :: libint2_build_eri
#if INCLUDE_ERI >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1), &
      BIND(C) :: libint2_build_eri1
#endif
#if INCLUDE_ERI >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2), &
      BIND(C) :: libint2_build_eri2
#endif
#endif

#ifdef INCLUDE_ERI2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri), &
      BIND(C) :: libint2_build_2eri
#if INCLUDE_ERI2 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1), &
      BIND(C) :: libint2_build_2eri1
#endif
#if INCLUDE_ERI2 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2), &
      BIND(C) :: libint2_build_2eri2
#endif
#endif

#ifdef INCLUDE_ERI3
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri), &
      BIND(C) :: libint2_build_3eri
#if INCLUDE_ERI3 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1), &
      BIND(C) :: libint2_build_3eri1
#endif
#if INCLUDE_ERI3 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2), &
      BIND(C) :: libint2_build_3eri2
#endif
#endif

   INTERFACE
      SUBROUTINE libint2_static_init() BIND(C)
      END SUBROUTINE

      SUBROUTINE libint2_static_cleanup() BIND(C)
      END SUBROUTINE

#ifdef INCLUDE_ERI
      SUBROUTINE libint2_init_eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri
      END FUNCTION

#if INCLUDE_ERI >= 1
      SUBROUTINE libint2_init_eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri1
      END FUNCTION
#endif

#if INCLUDE_ERI >= 2
      SUBROUTINE libint2_init_eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri2
      END FUNCTION
#endif
#endif

#ifdef INCLUDE_ERI2
      SUBROUTINE libint2_init_2eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri
      END FUNCTION

#if INCLUDE_ERI2 >= 1
      SUBROUTINE libint2_init_2eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri1
      END FUNCTION
#endif

#if INCLUDE_ERI2 >= 2
      SUBROUTINE libint2_init_2eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri2
      END FUNCTION
#endif
#endif

#ifdef INCLUDE_ERI3
      SUBROUTINE libint2_init_3eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri
      END FUNCTION

#if INCLUDE_ERI3 >= 1
      SUBROUTINE libint2_init_3eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri1
      END FUNCTION
#endif

#if INCLUDE_ERI3 >= 2
      SUBROUTINE libint2_init_3eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri2
      END FUNCTION
#endif
#endif
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE libint2_build(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE
   END INTERFACE

END MODULE
