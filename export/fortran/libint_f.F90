MODULE libint_f
   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_PTR, C_NULL_PTR, C_INT, C_FUNPTR, C_F_POINTER, C_F_PROCPOINTER, C_SIZE_T

#include <libint2/config.h>
#include <libint2/util/generated/libint2_params.h>
#include "fortran_incldefs.h"

   IMPLICIT NONE

#ifdef LIBINT2_MAX_AM
   INTEGER, PARAMETER :: libint2_max_am = LIBINT2_MAX_AM
#endif
#ifdef LIBINT2_MAX_AM_default
   INTEGER, PARAMETER :: libint2_max_am_default = LIBINT2_MAX_AM_default
#else
#  error "LIBINT2_MAX_AM_default is expected to be defined, libint2_params.h is misgenerated"
#endif
#ifdef LIBINT2_MAX_AM_default1
   INTEGER, PARAMETER :: libint2_max_am_default1 = LIBINT2_MAX_AM_default1
#else
   INTEGER, PARAMETER :: libint2_max_am_default1 = LIBINT2_MAX_AM_default
#endif
#ifdef LIBINT2_MAX_AM_default2
   INTEGER, PARAMETER :: libint2_max_am_default2 = LIBINT2_MAX_AM_default2
#else
   INTEGER, PARAMETER :: libint2_max_am_default2 = LIBINT2_MAX_AM_default
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
#ifdef LIBINT2_MAX_AM_3eri
   INTEGER, PARAMETER :: libint2_max_am_3eri = LIBINT2_MAX_AM_3eri
#endif
#ifdef LIBINT2_MAX_AM_3eri1
   INTEGER, PARAMETER :: libint2_max_am_3eri1 = LIBINT2_MAX_AM_3eri1
#endif
#ifdef LIBINT2_MAX_AM_3eri2
   INTEGER, PARAMETER :: libint2_max_am_3eri2 = LIBINT2_MAX_AM_3eri2
#endif
#ifdef LIBINT2_MAX_AM_2eri
   INTEGER, PARAMETER :: libint2_max_am_2eri = LIBINT2_MAX_AM_2eri
#endif
#ifdef LIBINT2_MAX_AM_2eri1
   INTEGER, PARAMETER :: libint2_max_am_2eri1 = LIBINT2_MAX_AM_2eri1
#endif
#ifdef LIBINT2_MAX_AM_2eri2
   INTEGER, PARAMETER :: libint2_max_am_2eri2 = LIBINT2_MAX_AM_2eri2
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
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri, 0:libint2_max_am_2eri), &
      BIND(C) :: libint2_build_2eri
#if INCLUDE_ERI2 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri1, 0:libint2_max_am_2eri1), &
      BIND(C) :: libint2_build_2eri1
#endif
#if INCLUDE_ERI2 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri2, 0:libint2_max_am_2eri2), &
      BIND(C) :: libint2_build_2eri2
#endif
#endif

#ifdef INCLUDE_ERI3
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default, 0:libint2_max_am_default, 0:libint2_max_am_3eri), &
      BIND(C) :: libint2_build_3eri
#if INCLUDE_ERI3 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default1, 0:libint2_max_am_default1, 0:libint2_max_am_3eri1), &
      BIND(C) :: libint2_build_3eri1
#endif
#if INCLUDE_ERI3 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default2, 0:libint2_max_am_default2, 0:libint2_max_am_3eri2), &
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
   SUBROUTINE compute_eri_f(contrdepth, deriv_order, am1, c1, alpha1, A, &
                            am2, c2, alpha2, B, &
                            am3, c3, alpha3, C, &
                            am4, c4, alpha4, D, F, erieval) bind(c)
      INTEGER, PARAMETER :: dp = C_DOUBLE
      INTEGER, PARAMETER :: i = C_INT
      INTEGER(KIND=i), INTENT(IN) :: contrdepth
      INTEGER(KIND=i), INTENT(IN) :: deriv_order
      INTEGER(KIND=i), INTENT(IN) :: am1, am2, am3, am4
      REAL(KIND=dp), DIMENSION(contrdepth), INTENT(IN) :: c1, c2, c3, c4
      REAL(KIND=dp), DIMENSION(contrdepth), INTENT(IN) :: alpha1, alpha2, alpha3, alpha4
      REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: A, B, C, D
      REAL(KIND=dp), DIMENSION(am1 + am2 + am3 + am4 + 1 + deriv_order, contrdepth**4), INTENT(IN) :: F
      TYPE(libint_t), DIMENSION(contrdepth**4), INTENT(OUT) :: erieval
      REAL(KIND=dp) :: gammap, AB2, CD2, PQ2, rhop, rhoq, gammaq, gammapq, K1, K2, pfac, C1234, &
                       alpha1_, alpha2_, alpha3_, alpha4_
      REAL(KIND=dp), DIMENSION(3) :: P, Q, QC, QD, PQ, PA, PB, W
      INTEGER(KIND=i) :: p1, p2, p3, p4, p1234, am, max_am, am_tot
      REAL(KIND=dp), PARAMETER :: pi = 4*atan(1.0_dp)
      PROCEDURE(libint2_build), POINTER :: build_eri

      am = am1 + am2 + am3 + am4
      max_am = MAXVAL([am1, am2, am3, am4])

      p1234 = 0
      DO p1 = 1, contrdepth
         DO p2 = 1, contrdepth
            DO p3 = 1, contrdepth
               DO p4 = 1, contrdepth

                  p1234 = p1234 + 1
                  alpha1_ = alpha1(p1)
                  alpha2_ = alpha2(p2)
                  alpha3_ = alpha3(p3)
                  alpha4_ = alpha4(p4)

                  gammap = alpha1_ + alpha2_
                  P = (alpha1_*A + alpha2_*B)/gammap
                  PA = P - A
                  PB = P - B
                  AB2 = SUM((A - B)*(A - B))

#if LIBINT2_DEFINED_PA_x
                  erieval(p1234)%PA_x(1) = PA(1)
#endif
#if LIBINT2_DEFINED_PA_y
                  erieval(p1234)%PA_y(1) = PA(2)
#endif
#if LIBINT2_DEFINED_PA_z
                  erieval(p1234)%PA_z(1) = PA(3)
#endif
#if LIBINT2_DEFINED_AB_x
                  erieval(p1234)%AB_x(1) = A(1) - B(1)
#endif
#if LIBINT2_DEFINED_AB_y
                  erieval(p1234)%AB_y(1) = A(2) - B(2)
#endif
#if LIBINT2_DEFINED_AB_z
                  erieval(p1234)%AB_z(1) = A(3) - B(3)
#endif
#if LIBINT2_DEFINED_PB_x
                  erieval(p1234)%PB_x(1) = PB(1)
#endif
#if LIBINT2_DEFINED_PB_y
                  erieval(p1234)%PB_y(1) = PB(2)
#endif
#if LIBINT2_DEFINED_PB_z
                  erieval(p1234)%PB_z(1) = PB(3)
#endif
#if LIBINT2_DEFINED_BA_x
                  erieval(p1234)%BA_x(1) = B(1) - A(1)
#endif
#if LIBINT2_DEFINED_BA_y
                  erieval(p1234)%BA_y(1) = B(2) - A(2)
#endif
#if LIBINT2_DEFINED_BA_z
                  erieval(p1234)%BA_z(1) = B(3) - A(3)
#endif
#if LIBINT2_DEFINED_oo2z
                  erieval(p1234)%oo2z(1) = 0.5_dp/gammap
#endif
                  gammaq = alpha3_ + alpha4_
                  gammapq = gammap*gammaq/(gammap + gammaq)
                  Q = (alpha3_*C + alpha4_*D)/gammaq
                  QC = Q - C
                  QD = Q - D
                  CD2 = SUM((C - D)*(C - D))

#if LIBINT2_DEFINED_QC_x
                  erieval(p1234)%QC_x(1) = QC(1)
#endif
#if LIBINT2_DEFINED_QC_y
                  erieval(p1234)%QC_y(1) = QC(2)
#endif
#if LIBINT2_DEFINED_QC_z
                  erieval(p1234)%QC_z(1) = QC(3)
#endif
#if LIBINT2_DEFINED_QD_x
                  erieval(p1234)%QD_x(1) = QD(1)
#endif
#if LIBINT2_DEFINED_QD_y
                  erieval(p1234)%QD_y(1) = QD(2)
#endif
#if LIBINT2_DEFINED_QD_z
                  erieval(p1234)%QD_z(1) = QD(3)
#endif
#if LIBINT2_DEFINED_CD_x
                  erieval(p1234)%CD_x(1) = C(1) - D(1)
#endif
#if LIBINT2_DEFINED_CD_y
                  erieval(p1234)%CD_y(1) = C(2) - D(2)
#endif
#if LIBINT2_DEFINED_CD_z
                  erieval(p1234)%CD_z(1) = C(3) - D(3)
#endif
#if LIBINT2_DEFINED_DC_x
                  erieval(p1234)%DC_x(1) = D(1) - C(1)
#endif
#if LIBINT2_DEFINED_DC_y
                  erieval(p1234)%DC_y(1) = D(2) - C(2)
#endif
#if LIBINT2_DEFINED_DC_z
                  erieval(p1234)%DC_z(1) = D(3) - C(3)
#endif
#if LIBINT2_DEFINED_oo2e
                  erieval(p1234)%oo2e(1) = 0.5_dp/gammaq
#endif
                  PQ = P - Q
                  PQ2 = SUM(PQ*PQ)
                  W = (gammap*P + gammaq*Q)/(gammap + gammaq)

#if LIBINT2_DEFINED_WP_x
                  erieval(p1234)%WP_x(1) = W(1) - P(1)
#endif
#if LIBINT2_DEFINED_WP_y
                  erieval(p1234)%WP_y(1) = W(2) - P(2)
#endif
#if LIBINT2_DEFINED_WP_z
                  erieval(p1234)%WP_z(1) = W(3) - P(3)
#endif
#if LIBINT2_DEFINED_WQ_x
                  erieval(p1234)%WQ_x(1) = W(1) - Q(1)
#endif
#if LIBINT2_DEFINED_WQ_y
                  erieval(p1234)%WQ_y(1) = W(2) - Q(2)
#endif
#if LIBINT2_DEFINED_WQ_z
                  erieval(p1234)%WQ_z(1) = W(3) - Q(3)
#endif
#if LIBINT2_DEFINED_oo2ze
                  erieval(p1234)%oo2ze(1) = 0.5_dp/(gammap + gammaq)
#endif

#if LIBINT2_DEFINED_roz
                  erieval(p1234)%roz(1) = gammapq/gammap
#endif
#if LIBINT2_DEFINED_roe
                  erieval(p1234)%roe(1) = gammapq/gammaq
#endif
                  IF (deriv_order > 0) THEN
#if LIBINT2_DEFINED_alpha1_rho_over_zeta2
                     erieval(p1234)%alpha1_rho_over_zeta2(1) = alpha1_*gammapq/(gammap*gammap)
#endif
#if LIBINT2_DEFINED_alpha2_rho_over_zeta2
                     erieval(p1234)%alpha2_rho_over_zeta2(1) = alpha2_*gammapq/(gammap*gammap)
#endif
#if LIBINT2_DEFINED_alpha3_rho_over_eta2
                     erieval(p1234)%alpha3_rho_over_eta2(1) = alpha3_*gammapq/(gammaq*gammaq)
#endif
#if LIBINT2_DEFINED_alpha4_rho_over_eta2
                     erieval(p1234)%alpha4_rho_over_eta2(1) = alpha4_*gammapq/(gammaq*gammaq)
#endif
#if LIBINT2_DEFINED_alpha1_over_zetapluseta
                     erieval(p1234)%alpha1_over_zetapluseta(1) = alpha1_/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha2_over_zetapluseta
                     erieval(p1234)%alpha2_over_zetapluseta(1) = alpha2_/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha3_over_zetapluseta
                     erieval(p1234)%alpha3_over_zetapluseta(1) = alpha3_/(gammap + gammaq)
#endif
#if LIBINT2_DEFINED_alpha4_over_zetapluseta
                     erieval(p1234)%alpha4_over_zetapluseta(1) = alpha4_/(gammap + gammaq)
#endif
                     rhop = alpha1_*alpha2_/gammap
#if LIBINT2_DEFINED_rho12_over_alpha1
                     erieval(p1234)%rho12_over_alpha1(1) = rhop/alpha1_
#endif
#if LIBINT2_DEFINED_rho12_over_alpha2
                     erieval(p1234)%rho12_over_alpha2(1) = rhop/alpha2_
#endif
                     rhoq = alpha3_*alpha4_/gammaq
#if LIBINT2_DEFINED_rho34_over_alpha3
                     erieval(p1234)%rho34_over_alpha3(1) = rhoq/alpha3_
#endif
#if LIBINT2_DEFINED_rho34_over_alpha4
                     erieval(p1234)%rho34_over_alpha4(1) = rhoq/alpha4_
#endif
#if LIBINT2_DEFINED_two_alpha0_bra
                     erieval(p1234)%two_alpha0_bra(1) = 2.0_dp*alpha1_
#endif
#if LIBINT2_DEFINED_two_alpha0_ket
                     erieval(p1234)%two_alpha0_ket(1) = 2.0_dp*alpha2_
#endif
#if LIBINT2_DEFINED_two_alpha1_ket
                     erieval(p1234)%two_alpha1_ket(1) = 2.0_dp*alpha4_
#endif
#if LIBINT2_DEFINED_two_alpha1_bra
                     erieval(p1234)%two_alpha1_bra(1) = 2.0_dp*alpha3_
#endif
                  ENDIF
                  K1 = exp(-alpha1_*alpha2_*AB2/gammap)
                  K2 = exp(-alpha3_*alpha4_*CD2/gammaq)
                  pfac = 2*pi**2.5_dp*K1*K2/(gammap*gammaq*SQRT(gammap + gammaq))
                  C1234 = c1(p1)*c2(p2)*c3(p3)*c4(p4)

                  am_tot = am + deriv_order

#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
                  IF (am_tot >= 0) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0(1) = &
                     C1234*pfac*F(1, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
                  IF (am_tot >= 1) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1(1) = &
                     C1234*pfac*F(2, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
                  IF (am_tot >= 2) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2(1) = &
                     C1234*pfac*F(3, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
                  IF (am_tot >= 3) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3(1) = &
                     C1234*pfac*F(4, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
                  IF (am_tot >= 4) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4(1) = &
                     C1234*pfac*F(5, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
                  IF (am_tot >= 5) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5(1) = &
                     C1234*pfac*F(6, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
                  IF (am_tot >= 6) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6(1) = &
                     C1234*pfac*F(7, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
                  IF (am_tot >= 7) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7(1) = &
                     C1234*pfac*F(8, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
                  IF (am_tot >= 8) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8(1) = &
                     C1234*pfac*F(9, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
                  IF (am_tot >= 9) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9(1) = &
                     C1234*pfac*F(10, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
                  IF (am_tot >= 10) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10(1) = &
                     C1234*pfac*F(11, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
                  IF (am_tot >= 11) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11(1) = &
                     C1234*pfac*F(12, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
                  IF (am_tot >= 12) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12(1) = &
                     C1234*pfac*F(13, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
                  IF (am_tot >= 13) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13(1) = &
                     C1234*pfac*F(14, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
                  IF (am_tot >= 14) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14(1) = &
                     C1234*pfac*F(15, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
                  IF (am_tot >= 15) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15(1) = &
                     C1234*pfac*F(16, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
                  IF (am_tot >= 16) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16(1) = &
                     C1234*pfac*F(17, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
                  IF (am_tot >= 17) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17(1) = &
                     C1234*pfac*F(18, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
                  IF (am_tot >= 18) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18(1) = &
                     C1234*pfac*F(19, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
                  IF (am_tot >= 19) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19(1) = &
                     C1234*pfac*F(20, p1234)
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
                  IF (am_tot >= 20) erieval(p1234)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20(1) = &
                     C1234*pfac*F(21, p1234)
#endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO

#if LIBINT2_CONTRACTED_INTS
      erieval(1)%contrdepth = p1234
#endif

      IF (deriv_order == 0) CALL C_F_PROCPOINTER(libint2_build_eri(am4, am3, am2, am1), build_eri)
#if INCLUDE_ERI >= 1
      IF (deriv_order == 1) CALL C_F_PROCPOINTER(libint2_build_eri1(am4, am3, am2, am1), build_eri)
#endif
#if INCLUDE_ERI >= 2
      IF (deriv_order == 2) CALL C_F_PROCPOINTER(libint2_build_eri2(am4, am3, am2, am1), build_eri)
#endif

      CALL build_eri(erieval)
   END SUBROUTINE
#endif

END MODULE
