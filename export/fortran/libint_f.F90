MODULE libint_f
   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_PTR, C_NULL_PTR, C_INT, C_FUNPTR, C_F_POINTER, C_F_PROCPOINTER, C_SIZE_T

#include <libint2/config.h>
#include <libint2/util/generated/libint2_params.h>
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

#include "libint2_types_f.h"

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

#ifdef INCLUDE_ERI
CONTAINS
   SUBROUTINE compute_eri(am1, alpha1, A, &
                          am2, alpha2, B, &
                          am3, alpha3, C, &
                          am4, alpha4, D, F, deriv_order) bind(c, name='compute_eri_f')
      INTEGER, INTENT(IN) :: am1, am2, am3, am4
      INTEGER, PARAMETER :: dp = C_DOUBLE
      REAL(KIND=dp), INTENT(IN) :: alpha1, alpha2, alpha3, alpha4
      REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: A, B, C, D
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN) :: F
      INTEGER, INTENT(IN) :: deriv_order

      REAL(KIND=dp), PARAMETER :: pi = 4*atan(1.0_dp)
      INTEGER, PARAMETER, DIMENSION(0:2) :: n_targets = [1, 12, 78]
      TYPE(libint_t), DIMENSION(1) :: erieval
      REAL(KIND=dp), DIMENSION(0:UBOUND(F, 1)) :: ssss
      INTEGER :: max_am, n1, n2, n3, n4, na, nb, nc, nd, ishell, i_target
      REAL(KIND=dp) :: gammap, AB2, CD2, PQ2, rhop, rhoq, gammaq, gammapq, K1, K2, pfac
      REAL(KIND=dp), DIMENSION(3) :: P, Q, QC, QD, PQ, PA, PB, W
      REAL(KIND=dp), DIMENSION(:), POINTER :: eri_shell_set
      PROCEDURE(libint2_build), POINTER :: build_eri

      max_am = MAXVAL([am1, am2, am3, am4])

      IF (deriv_order == 0) CALL libint2_init_eri(erieval, max_am, C_NULL_PTR)
#if INCLUDE_ERI >= 1
      IF (deriv_order == 1) CALL libint2_init_eri1(erieval, max_am, C_NULL_PTR)
#endif
#if INCLUDE_ERI >= 2
      IF (deriv_order == 2) CALL libint2_init_eri2(erieval, max_am, C_NULL_PTR)
#endif

#if LIBINT2_CONTRACTED_INTS
      erieval(1)%contrdepth = 1
#endif

      gammap = alpha1 + alpha2
      P = (alpha1*A + alpha2*B)/gammap
      PA = P - A
      PB = P - B
      AB2 = SUM((A - B)*(A - B))

#if LIBINT2_DEFINED_PA_x
      erieval(1)%PA_x(1) = PA(1)
#endif
#if LIBINT2_DEFINED_PA_y
      erieval(1)%PA_y(1) = PA(2)
#endif
#if LIBINT2_DEFINED_PA_z
      erieval(1)%PA_z(1) = PA(3)
#endif
#if LIBINT2_DEFINED_AB_x
      erieval(1)%AB_x(1) = A(1) - B(1)
#endif
#if LIBINT2_DEFINED_AB_y
      erieval(1)%AB_y(1) = A(2) - B(2)
#endif
#if LIBINT2_DEFINED_AB_z
      erieval(1)%AB_z(1) = A(3) - B(3)
#endif
#if LIBINT2_DEFINED_PB_x
      erieval(1)%PB_x(1) = PB(1)
#endif
#if LIBINT2_DEFINED_PB_y
      erieval(1)%PB_y(1) = PB(2)
#endif
#if LIBINT2_DEFINED_PB_z
      erieval(1)%PB_z(1) = PB(3)
#endif
#if LIBINT2_DEFINED_BA_x
      erieval(1)%BA_x(1) = B(1) - A(1)
#endif
#if LIBINT2_DEFINED_BA_y
      erieval(1)%BA_y(1) = B(2) - A(2)
#endif
#if LIBINT2_DEFINED_BA_z
      erieval(1)%BA_z(1) = B(3) - A(3)
#endif
#if LIBINT2_DEFINED_oo2z
      erieval(1)%oo2z(1) = 0.5_dp/gammap
#endif
      gammaq = alpha3 + alpha4
      gammapq = gammap*gammaq/(gammap + gammaq)
      Q = (alpha3*C + alpha4*D)/gammaq
      QC = Q - C
      QD = Q - D
      CD2 = SUM((C - D)*(C - D))
#if LIBINT2_DEFINED_QC_x
      erieval(1)%QC_x(1) = QC(1)
#endif
#if LIBINT2_DEFINED_QC_y
      erieval(1)%QC_y(1) = QC(2)
#endif
#if LIBINT2_DEFINED_QC_z
      erieval(1)%QC_z(1) = QC(3)
#endif
#if LIBINT2_DEFINED_QD_x
      erieval(1)%QD_x(1) = QD(1)
#endif
#if LIBINT2_DEFINED_QD_y
      erieval(1)%QD_y(1) = QD(2)
#endif
#if LIBINT2_DEFINED_QD_z
      erieval(1)%QD_z(1) = QD(3)
#endif
#if LIBINT2_DEFINED_CD_x
      erieval(1)%CD_x(1) = C(1) - D(1)
#endif
#if LIBINT2_DEFINED_CD_y
      erieval(1)%CD_y(1) = C(2) - D(2)
#endif
#if LIBINT2_DEFINED_CD_z
      erieval(1)%CD_z(1) = C(3) - D(3)
#endif
#if LIBINT2_DEFINED_DC_x
      erieval(1)%DC_x(1) = D(1) - C(1)
#endif
#if LIBINT2_DEFINED_DC_y
      erieval(1)%DC_y(1) = D(2) - C(2)
#endif
#if LIBINT2_DEFINED_DC_z
      erieval(1)%DC_z(1) = D(3) - C(3)
#endif
#if LIBINT2_DEFINED_oo2e
      erieval(1)%oo2e(1) = 0.5_dp/gammaq
#endif
      PQ = P - Q
      PQ2 = SUM(PQ*PQ)
      W = (gammap*P + gammaq*Q)/(gammap + gammaq)
#if LIBINT2_DEFINED_WP_x
      erieval(1)%WP_x(1) = W(1) - P(1)
#endif
#if LIBINT2_DEFINED_WP_y
      erieval(1)%WP_y(1) = W(2) - P(2)
#endif
#if LIBINT2_DEFINED_WP_z
      erieval(1)%WP_z(1) = W(3) - P(3)
#endif
#if LIBINT2_DEFINED_WQ_x
      erieval(1)%WQ_x(1) = W(1) - Q(1)
#endif
#if LIBINT2_DEFINED_WQ_y
      erieval(1)%WQ_y(1) = W(2) - Q(2)
#endif
#if LIBINT2_DEFINED_WQ_z
      erieval(1)%WQ_z(1) = W(3) - Q(3)
#endif
#if LIBINT2_DEFINED_oo2ze
      erieval(1)%oo2ze(1) = 0.5_dp/(gammap + gammaq)
#endif

#if LIBINT2_DEFINED_roz
      erieval(1)%roz(1) = gammapq/gammap
#endif
#if LIBINT2_DEFINED_roe
      erieval(1)%roe(1) = gammapq/gammaq
#endif
      IF (deriv_order > 0) THEN
#if LIBINT2_DEFINED_alpha1_rho_over_zeta2
         erieval(1)%alpha1_rho_over_zeta2(1) = alpha1*gammapq/(gammap*gammap)
#endif
#if LIBINT2_DEFINED_alpha2_rho_over_zeta2
         erieval(1)%alpha2_rho_over_zeta2(1) = alpha2*gammapq/(gammap*gammap)
#endif
#if LIBINT2_DEFINED_alpha3_rho_over_eta2
         erieval(1)%alpha3_rho_over_eta2(1) = alpha3*gammapq/(gammaq*gammaq)
#endif
#if LIBINT2_DEFINED_alpha4_rho_over_eta2
         erieval(1)%alpha4_rho_over_eta2(1) = alpha4*gammapq/(gammaq*gammaq)
#endif
#if LIBINT2_DEFINED_alpha1_over_zetapluseta
         erieval(1)%alpha1_over_zetapluseta(1) = alpha1/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha2_over_zetapluseta
         erieval(1)%alpha2_over_zetapluseta(1) = alpha2/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha3_over_zetapluseta
         erieval(1)%alpha3_over_zetapluseta(1) = alpha3/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha4_over_zetapluseta
         erieval(1)%alpha4_over_zetapluseta(1) = alpha4/(gammap + gammaq)
#endif
         rhop = alpha1*alpha2/gammap
#if LIBINT2_DEFINED_rho12_over_alpha1
         erieval(1)%rho12_over_alpha1(1) = rhop/alpha1
#endif
#if LIBINT2_DEFINED_rho12_over_alpha2
         erieval(1)%rho12_over_alpha2(1) = rhop/alpha2
#endif
         rhoq = alpha3*alpha4/gammaq
#if LIBINT2_DEFINED_rho34_over_alpha3
         erieval(1)%rho34_over_alpha3(1) = rhoq/alpha3
#endif
#if LIBINT2_DEFINED_rho34_over_alpha3
            erieval(1)%rho34_over_alpha3(1) = rhoq/alpha3
#endif
#if LIBINT2_DEFINED_rho34_over_alpha4
            erieval(1)%rho34_over_alpha4(1) = rhoq/alpha4
#endif
#if LIBINT2_DEFINED_two_alpha0_bra
         erieval(1)%two_alpha0_bra(1) = 2.0_dp*alpha1
#endif
#if LIBINT2_DEFINED_two_alpha0_ket
         erieval(1)%two_alpha0_ket(1) = 2.0_dp*alpha2
#endif
#if LIBINT2_DEFINED_two_alpha1_ket
         erieval(1)%two_alpha1_ket(1) = 2.0_dp*alpha4
#endif
#if LIBINT2_DEFINED_two_alpha1_bra
         erieval(1)%two_alpha1_bra(1) = 2.0_dp*alpha3
#endif
      ENDIF
      K1 = exp(-alpha1*alpha2*AB2/gammap)
      K2 = exp(-alpha3*alpha4*CD2/gammaq)
      pfac = 2*pi**2.5_dp*K1*K2/(gammap*gammaq*SQRT(gammap + gammaq))

      ssss = pfac*F

#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
      IF (UBOUND(ssss, 1) >= 0) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0(1) = ssss(0)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
      IF (UBOUND(ssss, 1) >= 1) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1(1) = ssss(1)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
      IF (UBOUND(ssss, 1) >= 2) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2(1) = ssss(2)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
      IF (UBOUND(ssss, 1) >= 3) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3(1) = ssss(3)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
      IF (UBOUND(ssss, 1) >= 4) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4(1) = ssss(4)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
      IF (UBOUND(ssss, 1) >= 5) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5(1) = ssss(5)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
      IF (UBOUND(ssss, 1) >= 6) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6(1) = ssss(6)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
      IF (UBOUND(ssss, 1) >= 7) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7(1) = ssss(7)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
      IF (UBOUND(ssss, 1) >= 8) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8(1) = ssss(8)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
      IF (UBOUND(ssss, 1) >= 9) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9(1) = ssss(9)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
      IF (UBOUND(ssss, 1) >= 10) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10(1) = ssss(10)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
      IF (UBOUND(ssss, 1) >= 11) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11(1) = ssss(11)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
      IF (UBOUND(ssss, 1) >= 12) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12(1) = ssss(12)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
      IF (UBOUND(ssss, 1) >= 13) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13(1) = ssss(13)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
      IF (UBOUND(ssss, 1) >= 14) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14(1) = ssss(14)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
      IF (UBOUND(ssss, 1) >= 15) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15(1) = ssss(15)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
      IF (UBOUND(ssss, 1) >= 16) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16(1) = ssss(16)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
      IF (UBOUND(ssss, 1) >= 17) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17(1) = ssss(17)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
      IF (UBOUND(ssss, 1) >= 17) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18(1) = ssss(18)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
      IF (UBOUND(ssss, 1) >= 17) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19(1) = ssss(19)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
      IF (UBOUND(ssss, 1) >= 17) erieval(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20(1) = ssss(20)
#endif

      IF (deriv_order == 0) CALL C_F_PROCPOINTER(libint2_build_eri(am4, am3, am2, am1), build_eri)
#if INCLUDE_ERI >= 1
      IF (deriv_order == 1) CALL C_F_PROCPOINTER(libint2_build_eri1(am4, am3, am2, am1), build_eri)
#endif
#if INCLUDE_ERI >= 2
      IF (deriv_order == 2) CALL C_F_PROCPOINTER(libint2_build_eri2(am4, am3, am2, am1), build_eri)
#endif

      CALL build_eri(erieval)

      n1 = (am1 + 1)*(am1 + 2)/2
      n2 = (am2 + 1)*(am2 + 2)/2
      n3 = (am3 + 1)*(am3 + 2)/2
      n4 = (am4 + 1)*(am4 + 2)/2

      WRITE (*, "(A14,I1)") "deriv order = ", deriv_order
      DO i_target = 1, n_targets(deriv_order)
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

      IF (deriv_order == 0) CALL libint2_cleanup_eri(erieval)
#if INCLUDE_ERI >= 1
      IF (deriv_order == 1) CALL libint2_cleanup_eri1(erieval)
#endif
#if INCLUDE_ERI >= 2
      IF (deriv_order == 2) CALL libint2_cleanup_eri2(erieval)
#endif

   END SUBROUTINE
#endif

END MODULE
