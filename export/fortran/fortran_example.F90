PROGRAM fortran_example

#include <libint2/config.h>
#include <libint2/util/generated/libint2_params.h>
#include "fortran_incldefs.h"

   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_F_POINTER, C_F_PROCPOINTER, C_NULL_PTR
   USE libint_f, ONLY: libint_t, libint2_static_init, libint2_static_cleanup, libint2_build, libint2_max_am_eri, &
      compute_eri_f

#ifdef INCLUDE_ERI
   USE libint_f, ONLY: libint2_init_eri, libint2_cleanup_eri
#if INCLUDE_ERI >= 1
   USE libint_f, ONLY: libint2_init_eri1, libint2_cleanup_eri1
#endif
#if INCLUDE_ERI >= 2
   USE libint_f, ONLY: libint2_init_eri2, libint2_cleanup_eri2
#endif
#endif

   IMPLICIT NONE

   INTEGER, PARAMETER :: dp = C_DOUBLE
   REAL(KIND=dp), DIMENSION(0:7) :: F
   REAL(KIND=dp), DIMENSION(3) :: A, B, C, D
   INTEGER :: max_am, am1, am2, am3, am4
   REAL(KIND=dp), DIMENSION(1) :: alpha1, alpha2, alpha3, alpha4
   REAL(KIND=dp), DIMENSION(1) :: c1, c2, c3, c4
   TYPE(libint_t), DIMENSION(1) :: erieval
   INTEGER :: deriv_order

   A = [0.1_dp, 1.3_dp, -1.5_dp]
   B = [2.0_dp, -0.5_dp, 1.1_dp]
   C = [-1.2_dp, 0.6_dp, -0.1_dp]
   D = [0.4_dp, 1.4_dp, 0.3_dp]

   am1 = MINVAL([1, libint2_max_am_eri])
   am2 = MINVAL([0, libint2_max_am_eri])
   am3 = MINVAL([2, libint2_max_am_eri])
   am4 = MINVAL([0, libint2_max_am_eri])
   max_am = MAXVAL([am1, am2, am3, am4])

   alpha1 = 0.1_dp
   alpha2 = 1.0_dp
   alpha3 = 0.5_dp
   alpha4 = 1.9_dp

   c1 = 1.0_dp
   c2 = 1.0_dp
   c3 = 1.0_dp
   c4 = 1.0_dp

   ! hard-coded values of Boys function, these must be calculated externally
   F(:) = [0.41608906765397796_dp, 0.044889937015574935_dp, &
           0.013706554295511562_dp, 0.0063780699489852013_dp, &
           0.39523364424416996_dp, 0.038762258204098128_dp, &
           0.010936175183774838_dp, 0.0047907366138884629_dp]

   CALL libint2_static_init()
#ifdef INCLUDE_ERI
   deriv_order = 0
   CALL libint2_init_eri(erieval, max_am, C_NULL_PTR)
   CALL compute_eri_f(1, deriv_order, am1, c1, alpha1, A, &
                    am2, c2, alpha2, B, &
                    am3, c3, alpha3, C, &
                    am4, c4, alpha4, D, &
                    F, erieval)
   CALL print_eri(am1, am2, am3, am4, deriv_order, erieval)
   CALL libint2_cleanup_eri(erieval)

#if INCLUDE_ERI >= 1
   deriv_order = 1
   CALL libint2_init_eri1(erieval, max_am, C_NULL_PTR)
   CALL compute_eri_f(1, deriv_order, am1, c1, alpha1, A, &
                    am2, c2, alpha2, B, &
                    am3, c3, alpha3, C, &
                    am4, c4, alpha4, D, &
                    F, erieval)
   CALL print_eri(am1, am2, am3, am4, deriv_order, erieval)
   CALL libint2_cleanup_eri1(erieval)
#endif
#if INCLUDE_ERI >= 2
   deriv_order = 2
   CALL libint2_init_eri2(erieval, max_am, C_NULL_PTR)
   CALL compute_eri_f(1, deriv_order, am1, c1, alpha1, A, &
                    am2, c2, alpha2, B, &
                    am3, c3, alpha3, C, &
                    am4, c4, alpha4, D, &
                    F, deriv_order, erieval)
   CALL print_eri(am1, am2, am3, am4, erieval)
   CALL libint2_cleanup_eri2(erieval)
#endif
#endif

   CALL libint2_static_cleanup()

CONTAINS
   SUBROUTINE print_eri(am1, am2, am3, am4, deriv_order, erieval)
      INTEGER, INTENT(IN) :: am1, am2, am3, am4, deriv_order
      TYPE(libint_t), DIMENSION(*), INTENT(IN) :: erieval
      REAL(KIND=dp), DIMENSION(:), POINTER :: eri_shell_set
      INTEGER :: n1, n2, n3, n4, i_target, na, nb, nc, nd, ishell
      INTEGER, PARAMETER, DIMENSION(3) :: n_targets = [1, 12, 78]

      n1 = (am1 + 1)*(am1 + 2)/2
      n2 = (am2 + 1)*(am2 + 2)/2
      n3 = (am3 + 1)*(am3 + 2)/2
      n4 = (am4 + 1)*(am4 + 2)/2

      WRITE (*, "(A14,I1)") "deriv order = ", deriv_order
      DO i_target = 1, n_targets(deriv_order + 1)
         WRITE (*, "(A5,1X,I2)") "Shell-set #", i_target
         CALL C_F_POINTER(erieval(1)%targets(i_target), eri_shell_set, SHAPE=[n1*n2*n3*n4])

         ishell = 0
         DO na = 1, n1
            DO nb = 1, n2
               DO nc = 1, n3
                  DO nd = 1, n4
                     ishell = ishell + 1
                     WRITE (*, "(A5, I4, A11, F21.17)") &
                        "Elem ", ishell, ", (ab|cd) = ", eri_shell_set(ishell)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE

END PROGRAM
