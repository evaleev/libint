#include "../tests/unit/catch.hpp"
#include "../tests/eri/prep_libint2.h"
#include <libint2/config.h>
#include <libint2.h>
#include <libint2/util/generated/libint2_params.h>
#include <libint2/util/generated/libint2_types.h>
#include <libint2/util/memory.h>
#include <libint2/boys.h>
#include <iomanip>


#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif


#ifdef INCLUDE_ERI
/// compute_eri from Fortran
void compute_eri_f(int &contrdepth, int &deriv_order, int &am1, double *c1, double *alpha1, double *A,
                   int &am2, double *c2, double *alpha2, double *B,
                   int &am3, double *c3, double *alpha3, double *C,
                   int &am4, double *c4, double *alpha4, double *D,
                   double *F, Libint_t *erieval);

/// Independent C++ implementation of compute_eri, in order to check Fortran result.
void compute_eri_c(int &contrdepth, int &deriv_order, int &am1, double *c1, double *alpha1, double *A,
                   int &am2, double *c2, double *alpha2, double *B,
                   int &am3, double *c3, double *alpha3, double *C,
                   int &am4, double *c4, double *alpha4, double *D,
                   double *F, Libint_t *erieval)
{

  const unsigned int am = am1+am2+am3+am4;

  int p0123 = 0;
  for(int p0=0; p0<contrdepth; ++p0){
    for(int p1=0; p1<contrdepth; ++p1){
      for(int p2=0; p2<contrdepth; ++p2){
        for(int p3=0; p3<contrdepth; ++p3){
          double alpha1_ = alpha1[p0];
          double alpha2_ = alpha2[p1];
          double alpha3_ = alpha3[p2];
          double alpha4_ = alpha4[p3];

          double gammap = alpha1_ + alpha2_;
          double Px = (alpha1_*A[0] + alpha2_*B[0])/gammap;
          double Py = (alpha1_*A[1] + alpha2_*B[1])/gammap;
          double Pz = (alpha1_*A[2] + alpha2_*B[2])/gammap;
          double PAx = Px - A[0];
          double PAy = Py - A[1];
          double PAz = Pz - A[2];
          double PBx = Px - B[0];
          double PBy = Py - B[1];
          double PBz = Pz - B[2];
          double AB2 = (A[0]-B[0])*(A[0]-B[0])
                     + (A[1]-B[1])*(A[1]-B[1])
                     + (A[2]-B[2])*(A[2]-B[2]);

#if LIBINT2_DEFINED_PA_x
          erieval[p0123].PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED_PA_y
          erieval[p0123].PA_y[0] = PAy;
#endif

#if LIBINT2_DEFINED_PA_z
          erieval[p0123].PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED_AB_x
          erieval[p0123].AB_x[0] = A[0] - B[0];
#endif
#if LIBINT2_DEFINED_AB_y
          erieval[p0123].AB_y[0] = A[1] - B[1];
#endif
#if LIBINT2_DEFINED_AB_z
          erieval[p0123].AB_z[0] = A[2] - B[2];
#endif
#if LIBINT2_DEFINED_PB_x
          erieval[p0123].PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED_PB_y
          erieval[p0123].PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED_PB_z
          erieval[p0123].PB_z[0] = PBz;
#endif
#if LIBINT2_DEFINED_BA_x
          erieval[p0123].BA_x[0] = B[0] - A[0];
#endif
#if LIBINT2_DEFINED_BA_y
          erieval[p0123].BA_y[0] = B[1] - A[1];
#endif
#if LIBINT2_DEFINED_BA_z
          erieval[p0123].BA_z[0] = B[2] - A[2];
#endif
#if LIBINT2_DEFINED_oo2z
          erieval[p0123].oo2z[0] = 0.5/gammap;
#endif
          double gammaq = alpha3_ + alpha4_;
          double gammapq = gammap*gammaq/(gammap+gammaq);
          double Qx = (alpha3_*C[0] + alpha4_*D[0])/gammaq;
          double Qy = (alpha3_*C[1] + alpha4_*D[1])/gammaq;
          double Qz = (alpha3_*C[2] + alpha4_*D[2])/gammaq;
          double QCx = Qx - C[0];
          double QCy = Qy - C[1];
          double QCz = Qz - C[2];
          double QDx = Qx - D[0];
          double QDy = Qy - D[1];
          double QDz = Qz - D[2];
          double CD2 = (C[0]-D[0])*(C[0]-D[0])
                     + (C[1]-D[1])*(C[1]-D[1])
                     + (C[2]-D[2])*(C[2]-D[2]);
#if LIBINT2_DEFINED_QC_x
          erieval[p0123].QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED_QC_y
          erieval[p0123].QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED_QC_z
          erieval[p0123].QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED_QD_x
          erieval[p0123].QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED_QD_y
          erieval[p0123].QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED_QD_z
          erieval[p0123].QD_z[0] = QDz;
#endif
#if LIBINT2_DEFINED_CD_x
          erieval[p0123].CD_x[0] = C[0] - D[0];
#endif
#if LIBINT2_DEFINED_CD_y
          erieval[p0123].CD_y[0] = C[1] - D[1];
#endif
#if LIBINT2_DEFINED_CD_z
          erieval[p0123].CD_z[0] = C[2] - D[2];
#endif
#if LIBINT2_DEFINED_DC_x
          erieval[p0123].DC_x[0] = D[0] - C[0];
#endif
#if LIBINT2_DEFINED_DC_y
          erieval[p0123].DC_y[0] = D[1] - C[1];
#endif
#if LIBINT2_DEFINED_DC_z
          erieval[p0123].DC_z[0] = D[2] - C[2];
#endif
#if LIBINT2_DEFINED_oo2e
          erieval[p0123].oo2e[0] = 0.5/gammaq;
#endif

          double PQx = Px - Qx;
          double PQy = Py - Qy;
          double PQz = Pz - Qz;
          double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
          double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
          double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
          double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

#if LIBINT2_DEFINED_WP_x
          erieval[p0123].WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED_WP_y
          erieval[p0123].WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED_WP_z
          erieval[p0123].WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED_WQ_x
          erieval[p0123].WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED_WQ_y
          erieval[p0123].WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED_WQ_z
          erieval[p0123].WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED_oo2ze
          erieval[p0123].oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED_roz
          erieval[p0123].roz[0] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED_roe
          erieval[p0123].roe[0] = gammapq/gammaq;
#endif
          if (deriv_order > 0){
#if LIBINT2_DEFINED_alpha1_rho_over_zeta2
            erieval[p0123].alpha1_rho_over_zeta2[0] = alpha1_*gammapq/(gammap*gammap);
#endif
#if LIBINT2_DEFINED_alpha2_rho_over_zeta2
            erieval[p0123].alpha2_rho_over_zeta2[0] = alpha2_*gammapq/(gammap*gammap);
#endif
#if LIBINT2_DEFINED_alpha3_rho_over_eta2
            erieval[p0123].alpha3_rho_over_eta2[0] = alpha3_*gammapq/(gammaq*gammaq);
#endif
#if LIBINT2_DEFINED_alpha4_rho_over_eta2
            erieval[p0123].alpha4_rho_over_eta2[0] = alpha4_*gammapq/(gammaq*gammaq);
#endif
#if LIBINT2_DEFINED_alpha1_over_zetapluseta
            erieval[p0123].alpha1_over_zetapluseta[0] = alpha1_/(gammap + gammaq);
#endif
#if LIBINT2_DEFINED_alpha2_over_zetapluseta
            erieval[p0123].alpha2_over_zetapluseta[0] = alpha2_/(gammap + gammaq);
#endif
#if LIBINT2_DEFINED_alpha3_over_zetapluseta
            erieval[p0123].alpha3_over_zetapluseta[0] = alpha3_/(gammap + gammaq);
#endif
#if LIBINT2_DEFINED_alpha4_over_zetapluseta
            erieval[p0123].alpha4_over_zetapluseta[0] = alpha4_/(gammap + gammaq);
#endif
            double rhop = alpha1_*alpha2_/gammap;
#if LIBINT2_DEFINED_rho12_over_alpha1
            erieval[p0123].rho12_over_alpha1[0] = rhop/alpha1_;
#endif
#if LIBINT2_DEFINED_rho12_over_alpha2
            erieval[p0123].rho12_over_alpha2[0] = rhop/alpha2_;
#endif
            double rhoq = alpha3_*alpha4_/gammaq;
#if LIBINT2_DEFINED_rho34_over_alpha3
            erieval[p0123].rho34_over_alpha3[0] = rhoq/alpha3_;
#endif
#if LIBINT2_DEFINED_rho34_over_alpha3
            erieval[p0123].rho34_over_alpha3[0] = rhoq/alpha3_;
#endif
#if LIBINT2_DEFINED_rho34_over_alpha4
            erieval[p0123].rho34_over_alpha4[0] = rhoq/alpha4_;
#endif
#if LIBINT2_DEFINED_two_alpha0_bra
            erieval[p0123].two_alpha0_bra[0] = 2.0*alpha1_;
#endif
#if LIBINT2_DEFINED_two_alpha0_ket
            erieval[p0123].two_alpha0_ket[0] = 2.0*alpha2_;
#endif
#if LIBINT2_DEFINED_two_alpha1_ket
            erieval[p0123].two_alpha1_ket[0] = 2.0*alpha4_;
#endif
#if LIBINT2_DEFINED_two_alpha1_bra
            erieval[p0123].two_alpha1_bra[0] = 2.0*alpha3_;
#endif
          }

          double K1 = exp(-alpha1_*alpha2_*AB2/gammap);
          double K2 = exp(-alpha3_*alpha4_*CD2/gammaq);
          double pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));
          double c1234 = c1[p0]*c2[p1]*c3[p2]*c4[p3];

          const int am_tot = am+deriv_order;

#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0
          if(am_tot >= 0) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0[0] = c1234*pfac*F[p0123*(am_tot + 1) + 0];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1
          if(am_tot >= 1) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[0] = c1234*pfac*F[p0123*(am_tot + 1) + 1];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2
          if(am_tot >= 2) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[0] = c1234*pfac*F[p0123*(am_tot + 1) + 2];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3
          if(am_tot >= 3) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3[0] = c1234*pfac*F[p0123*(am_tot + 1) + 3];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4
          if(am_tot >= 4) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4[0] = c1234*pfac*F[p0123*(am_tot + 1) + 4];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5
          if(am_tot >= 5) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5[0] = c1234*pfac*F[p0123*(am_tot + 1) + 5];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6
          if(am_tot >= 6) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6[0] = c1234*pfac*F[p0123*(am_tot + 1) + 6];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7
          if(am_tot >= 7) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7[0] = c1234*pfac*F[p0123*(am_tot + 1) + 7];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8
          if(am_tot >= 8) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8[0] = c1234*pfac*F[p0123*(am_tot + 1) + 8];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9
          if(am_tot >= 9) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9[0] = c1234*pfac*F[p0123*(am_tot + 1) + 9];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10
          if(am_tot >= 10) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10[0] = c1234*pfac*F[p0123*(am_tot + 1) + 10];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11
          if(am_tot >= 11) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11[0] = c1234*pfac*F[p0123*(am_tot + 1) + 11];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12
          if(am_tot >= 12) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12[0] = c1234*pfac*F[p0123*(am_tot + 1) + 12];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13
          if(am_tot >= 13) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13[0] = c1234*pfac*F[p0123*(am_tot + 1) + 13];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14
          if(am_tot >= 14) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14[0] = c1234*pfac*F[p0123*(am_tot + 1) + 14];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15
          if(am_tot >= 15) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15[0] = c1234*pfac*F[p0123*(am_tot + 1) + 15];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16
          if(am_tot >= 16) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16[0] = c1234*pfac*F[p0123*(am_tot + 1) + 16];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17
          if(am_tot >= 17) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17[0] = c1234*pfac*F[p0123*(am_tot + 1) + 17];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
          if(am_tot >= 18) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18[0] = c1234*pfac*F[p0123*(am_tot + 1) + 18];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19
          if(am_tot >= 19) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19[0] = c1234*pfac*F[p0123*(am_tot + 1) + 19];
#endif
#ifdef LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20
          if(am_tot >= 20) erieval[p0123]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20[0] = c1234*pfac*F[p0123*(am_tot + 1) + 20];
#endif

          ++p0123;
        }
      }
    }
  }

#if LIBINT2_CONTRACTED_INTS
  erieval[0].contrdepth = p0123;
#endif


  if (deriv_order == 0) libint2_build_eri[am1][am2][am3][am4](erieval);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) libint2_build_eri1[am1][am2][am3][am4](erieval);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) libint2_build_eri2[am1][am2][am3][am4](erieval);
#endif
}

bool test_eri_c_f(int &contrdepth, int &am1, double *c1, double *alpha1, double *A,
                     int &am2, double *c2, double *alpha2, double *B,
                     int &am3, double *c3, double *alpha3, double *C,
                     int &am4, double *c4, double *alpha4, double *D,
                     int &deriv_order, const double &threshold)
{
  const int am = am1 + am2 + am3 + am4;
  const int contrdepth4 = pow(contrdepth, 4);

  bool success = true;

  Libint_t* erieval_c = libint2::malloc<Libint_t>(contrdepth4);
  Libint_t* erieval_f = libint2::malloc<Libint_t>(contrdepth4);

  const int n_targets[] = {1, 12, 78};
  const unsigned int max_am = std::max(std::max(am1,am2),std::max(am3,am4));

  if (deriv_order == 0){
    libint2_init_eri(erieval_c, max_am, 0);
    libint2_init_eri(erieval_f, max_am, 0);
  }
#if INCLUDE_ERI >= 1
  if (deriv_order == 1){
    libint2_init_eri1(erieval_c, max_am, 0);
    libint2_init_eri1(erieval_f, max_am, 0);
  }
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2){
    libint2_init_eri2(erieval_c, max_am, 0);
    libint2_init_eri2(erieval_f, max_am, 0);
  }
#endif


  // Compute Boys function values

  const int am_tot = am + deriv_order;
  double F[contrdepth4*(am_tot+1)];
  libint2::FmEval_Chebyshev7<double> fmeval(am_tot);

  int p0123 = 0;
  for(int p0=0; p0<contrdepth; ++p0){
    for(int p1=0; p1<contrdepth; ++p1){
      for(int p2=0; p2<contrdepth; ++p2){
        for(int p3=0; p3<contrdepth; ++p3){
          double alpha1_ = alpha1[p0];
          double alpha2_ = alpha2[p1];
          double alpha3_ = alpha3[p2];
          double alpha4_ = alpha4[p3];

          double gammap = alpha1_ + alpha2_;
          double gammaq = alpha3_ + alpha4_;
          double gammapq = gammap*gammaq/(gammap+gammaq);

          double Qx = (alpha3_*C[0] + alpha4_*D[0])/gammaq;
          double Qy = (alpha3_*C[1] + alpha4_*D[1])/gammaq;
          double Qz = (alpha3_*C[2] + alpha4_*D[2])/gammaq;

          double Px = (alpha1_*A[0] + alpha2_*B[0])/gammap;
          double Py = (alpha1_*A[1] + alpha2_*B[1])/gammap;
          double Pz = (alpha1_*A[2] + alpha2_*B[2])/gammap;

          double PQx = Px - Qx;
          double PQy = Py - Qy;
          double PQz = Pz - Qz;

          double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
          fmeval.eval(&F[p0123*(am_tot+1)], PQ2 * gammapq, am_tot);

          ++p0123;
        }
      }
    }
  }

  // Compute integrals using C++ implementation
  compute_eri_c(contrdepth, deriv_order, am1, c1, alpha1, A,
                am2, c2, alpha2, B,
                am3, c3, alpha3, C,
                am4, c4, alpha4, D, F, erieval_c);

  // Compute integrals using Fortran implementation
  compute_eri_f(contrdepth, deriv_order, am1, c1, alpha1, A,
                am2, c2, alpha2, B,
                am3, c3, alpha3, C,
                am4, c4, alpha4, D, F, erieval_f);

  // Compare results
  for(int i_target=0; i_target<n_targets[deriv_order]; i_target++) {

    const double* eri_shell_set_c = erieval_c[0].targets[i_target];
    const double* eri_shell_set_f = erieval_f[0].targets[i_target];

    const unsigned int n1 = (am1 + 1) * (am1 + 2)/2;
    const unsigned int n2 = (am2 + 1) * (am2 + 2)/2;
    const unsigned int n3 = (am3 + 1) * (am3 + 2)/2;
    const unsigned int n4 = (am4 + 1) * (am4 + 2)/2;
    unsigned int nel = 0;
    for(int a=0; a<n1; a++) {
      for(int b=0; b<n2; b++) {
        for(int c=0; c<n3; c++) {
          for(int d=0; d<n4; d++, nel++) {
            const double abs_error = abs(*eri_shell_set_c - *eri_shell_set_f);
            if(abs_error > threshold) {
              std::cout << std::setprecision(17) << "Elem " << nel << " di= " << deriv_order <<
                ", : C = " << *eri_shell_set_c
                << ", Fortran = " << *eri_shell_set_f
                << ", abs_error = " << abs_error << std::endl;
              success = false;
            }
            ++eri_shell_set_c; ++eri_shell_set_f;
          }
        }
      }
    }
  }

  if (deriv_order == 0) {
     libint2_cleanup_eri(erieval_c);
     libint2_cleanup_eri(erieval_f);
   }
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) {
     libint2_cleanup_eri1(erieval_c);
     libint2_cleanup_eri1(erieval_f);
  }
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) {
     libint2_cleanup_eri2(erieval_c);
     libint2_cleanup_eri2(erieval_f);
  }
#endif

  free(erieval_c);
  free(erieval_f);

  return success;

}
#endif // INCLUDE_ERI

#if LIBINT2_SUPPORT_ERI

TEST_CASE("Fortran ERI", "[eri]") {

  const double threshold = 1.0E-16;

  int am1 = std::min(1, LIBINT2_MAX_AM_eri);
  int am2 = std::min(0, LIBINT2_MAX_AM_eri);
  int am3 = std::min(2, LIBINT2_MAX_AM_eri);
  int am4 = std::min(0, LIBINT2_MAX_AM_eri);

  unsigned int am[] = {(unsigned int)am1, (unsigned int)am2, (unsigned int)am3, (unsigned int)am4};

  int contrdepth = 2;
  unsigned int veclen = 1;

  unsigned int contrdepth_u = contrdepth;
  RandomShellSet<4u> rsqset(am, veclen, contrdepth_u);

  double* A = &(rsqset.R[0][0]);
  double* B = &(rsqset.R[1][0]);
  double* C = &(rsqset.R[2][0]);
  double* D = &(rsqset.R[3][0]);

  double alpha1[] = {rsqset.exp[0][0][0], rsqset.exp[0][0][1]};
  double alpha2[] = {rsqset.exp[1][0][0], rsqset.exp[1][0][1]};
  double alpha3[] = {rsqset.exp[2][0][0], rsqset.exp[2][0][1]};
  double alpha4[] = {rsqset.exp[3][0][0], rsqset.exp[3][0][1]};

  double c1[] = {rsqset.coef[0][0][0], rsqset.coef[0][0][1]};
  double c2[] = {rsqset.coef[1][0][0], rsqset.coef[1][0][1]};
  double c3[] = {rsqset.coef[2][0][0], rsqset.coef[2][0][1]};
  double c4[] = {rsqset.coef[3][0][0], rsqset.coef[3][0][1]};

#ifdef INCLUDE_ERI
   int deriv_order = 0;
   REQUIRE(test_eri_c_f(contrdepth, am1, c1, alpha1, A,
                           am2, c2, alpha2, B,
                           am3, c3, alpha3, C,
                           am4, c4, alpha4, D,
                           deriv_order, threshold));

#if INCLUDE_ERI >= 1
   deriv_order = 1;
   REQUIRE(test_eri_c_f(contrdepth, am1, c1, alpha1, A,
                           am2, c2, alpha2, B,
                           am3, c3, alpha3, C,
                           am4, c4, alpha4, D,
                           deriv_order, threshold));
#endif

#if INCLUDE_ERI >= 2
  deriv_order = 2;
  REQUIRE(test_eri_c_f(contrdepth, am1, c1, alpha1, A,
                          am2, c2, alpha2, B,
                          am3, c3, alpha3, C,
                          am4, c4, alpha4, D,
                          deriv_order, threshold));
#endif
#endif // INCLUDE_ERI
}

#endif // LIBINT2_SUPPORT_ERI
