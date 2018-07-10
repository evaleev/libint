PROGRAM fortran_example

#include <libint2/config.h>
#include <libint2/util/generated/libint2_params.h>
#include "fortran_incldefs.h"

   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_F_POINTER, C_F_PROCPOINTER, C_NULL_PTR
   USE libint_f, ONLY: libint_t, libint2_static_init, libint2_static_cleanup, libint2_build, libint2_max_am_eri, &
                       compute_eri

#ifdef INCLUDE_ERI
   USE libint_f, ONLY: libint2_init_eri, libint2_cleanup_eri, libint2_build_eri
#if INCLUDE_ERI >= 1
   USE libint_f, ONLY: libint2_init_eri1, libint2_cleanup_eri1, libint2_build_eri1
#endif
#if INCLUDE_ERI >= 2
   USE libint_f, ONLY: libint2_init_eri2, libint2_cleanup_eri2, libint2_build_eri2
#endif
#endif

   IMPLICIT NONE

   INTEGER, PARAMETER :: dp = C_DOUBLE
   REAL(KIND=dp), DIMENSION(0:7) :: F
   REAL(KIND=dp), DIMENSION(3) :: A, B, C, D
   INTEGER :: n1, n2, n3, n4
   REAL(KIND=dp) :: alpha1, alpha2, alpha3, alpha4

   A = [0.1_dp, 1.3_dp, -1.5_dp]
   B = [2.0_dp, -0.5_dp, 1.1_dp]
   C = [-1.2_dp, 0.6_dp, -0.1_dp]
   D = [0.4_dp, 1.4_dp, 0.3_dp]

   n1 = MINVAL([1, libint2_max_am_eri])
   n2 = MINVAL([0, libint2_max_am_eri])
   n3 = MINVAL([2, libint2_max_am_eri])
   n4 = MINVAL([0, libint2_max_am_eri])

   alpha1 = 0.1_dp
   alpha2 = 1.0_dp
   alpha3 = 0.5_dp
   alpha4 = 1.9_dp

   ! hard-coded values of Boys function for this input
   F = [0.41608906765397796_dp, 0.044889937015574935_dp, &
        0.013706554295511562_dp, 0.0063780699489852013_dp, &
        0.39523364424416996_dp, 0.038762258204098128_dp, &
        0.010936175183774838_dp, 0.0047907366138884629_dp]

   CALL libint2_static_init()
#ifdef INCLUDE_ERI
   CALL compute_eri(n1, alpha1, A, &
                    n2, alpha2, B, &
                    n3, alpha3, C, &
                    n4, alpha4, D, &
                    F, deriv_order=0)

#if INCLUDE_ERI >= 1
   CALL compute_eri(n1, alpha1, A, &
                    n2, alpha2, B, &
                    n3, alpha3, C, &
                    n4, alpha4, D, &
                    F, deriv_order=1)
#endif
#if INCLUDE_ERI >= 2
   CALL compute_eri(n1, alpha1, A, &
                    n2, alpha2, B, &
                    n3, alpha3, C, &
                    n4, alpha4, D, &
                    F, deriv_order=2)
#endif
#endif

   CALL libint2_static_cleanup()

END PROGRAM
