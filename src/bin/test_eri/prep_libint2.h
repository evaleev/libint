/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <libint2.h>
#include <libint2/boys.h>
#include <time.h>

#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

extern libint2::FmEval_Chebyshev7<double> fmeval_chebyshev;
extern libint2::FmEval_Taylor<double, 6> fmeval_taylor;

typedef unsigned int uint;
template <unsigned int N>
struct RandomShellSet {
 public:
  RandomShellSet(uint* am, uint veclen, uint contrdepth,
                 unsigned int seed = std::random_device{}()) {
    std::copy(am, am + N, l);

    std::mt19937 rng(seed);  // produces randomness out of thin air
    std::uniform_real_distribution<> rdist(0.7,
                                           1.3);  // distribution of coordinates
    std::uniform_real_distribution<> ecdist(
        0.1, 3.0);  // distribution of exponents/coefficients
    // glues source randomness with mapping
    auto rdie = [&rng, &rdist]() -> double { return rdist(rng); };
    auto ecdie = [&rng, &ecdist]() -> double { return ecdist(rng); };

    for (uint c = 0; c < N; ++c) {
      R[c].resize(3);
      std::generate(R[c].begin(), R[c].end(), rdie);

      exp[c].resize(veclen);
      coef[c].resize(veclen);
      for (uint v = 0; v < veclen; ++v) {
        exp[c][v].resize(contrdepth);
        std::generate(exp[c][v].begin(), exp[c][v].end(), ecdie);
        coef[c][v].resize(contrdepth);
        std::generate(coef[c][v].begin(), coef[c][v].end(), ecdie);
      }
    }
  }

  uint l[N];                                  // angular momenta
  std::vector<double> R[N];                   // origins
  std::vector<std::vector<double> > exp[N];   // exponents
  std::vector<std::vector<double> > coef[N];  // coefficients
};

template <typename LibintEval>
void prep_libint2(LibintEval* erievals, const RandomShellSet<4u>& rsqset,
                  int norm_flag, int deriv_order = 0) {
  const uint veclen = rsqset.exp[0].size();
  const uint contrdepth = rsqset.exp[0][0].size();

  const double* A = &(rsqset.R[0][0]);
  const double* B = &(rsqset.R[1][0]);
  const double* C = &(rsqset.R[2][0]);
  const double* D = &(rsqset.R[3][0]);

  const uint* am = rsqset.l;
  const unsigned int am01 = am[0] + am[1];
  const unsigned int am23 = am[2] + am[3];
  const unsigned int amtot = am01 + am23 + deriv_order;
  // static seems to be important for performance on OS X 10.7 with Apple
  // clang 3.1 (318.0.58) 8/24/2012 also seems to affect other compilers
  // (macports g++ 4.7.1) on OS X 10.7 and 10.8
  static double F[LIBINT_MAX_AM * 4 + 6];

  uint p0123 = 0;
  for (uint p0 = 0; p0 < contrdepth; p0++) {
    for (uint p1 = 0; p1 < contrdepth; p1++) {
      for (uint p2 = 0; p2 < contrdepth; p2++) {
        for (uint p3 = 0; p3 < contrdepth; p3++, p0123++) {
          LibintEval* erieval = &erievals[p0123];
          erieval->veclen = veclen;
#if LIBINT2_FLOP_COUNT
          erieval->nflops = erievals[0].nflops;
#endif

          for (uint v = 0; v < veclen; v++) {
            const double alpha0 = rsqset.exp[0][v][p0];
            const double alpha1 = rsqset.exp[1][v][p1];
            const double alpha2 = rsqset.exp[2][v][p2];
            const double alpha3 = rsqset.exp[3][v][p3];

            const double c0 = rsqset.coef[0][v][p0];
            const double c1 = rsqset.coef[1][v][p1];
            const double c2 = rsqset.coef[2][v][p2];
            const double c3 = rsqset.coef[3][v][p3];

            const double gammap = alpha0 + alpha1;
            const double oogammap = 1.0 / gammap;
            const double rhop = alpha0 * alpha1 * oogammap;
            const double Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
            const double Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
            const double Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
            const double AB_x = A[0] - B[0];
            const double AB_y = A[1] - B[1];
            const double AB_z = A[2] - B[2];
            const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

            const double PAx = Px - A[0];
            const double PAy = Py - A[1];
            const double PAz = Pz - A[2];
            const double PBx = Px - B[0];
            const double PBy = Py - B[1];
            const double PBz = Pz - B[2];

#if LIBINT2_DEFINED(eri, PA_x)
            erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri, PA_y)
            erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri, PA_z)
            erieval->PA_z[v] = PAz;
#endif
#if LIBINT2_DEFINED(eri, PB_x)
            erieval->PB_x[v] = PBx;
#endif
#if LIBINT2_DEFINED(eri, PB_y)
            erieval->PB_y[v] = PBy;
#endif
#if LIBINT2_DEFINED(eri, PB_z)
            erieval->PB_z[v] = PBz;
#endif

#if LIBINT2_DEFINED(eri, AB_x)
            erieval->AB_x[v] = AB_x;
#endif
#if LIBINT2_DEFINED(eri, AB_y)
            erieval->AB_y[v] = AB_y;
#endif
#if LIBINT2_DEFINED(eri, AB_z)
            erieval->AB_z[v] = AB_z;
#endif
#if LIBINT2_DEFINED(eri, BA_x)
            erieval->BA_x[v] = -AB_x;
#endif
#if LIBINT2_DEFINED(eri, BA_y)
            erieval->BA_y[v] = -AB_y;
#endif
#if LIBINT2_DEFINED(eri, BA_z)
            erieval->BA_z[v] = -AB_z;
#endif
#if LIBINT2_DEFINED(eri, oo2z)
            erieval->oo2z[v] = 0.5 * oogammap;
#endif

            const double gammaq = alpha2 + alpha3;
            const double oogammaq = 1.0 / gammaq;
            const double rhoq = alpha2 * alpha3 * oogammaq;
            const double one_o_gammap_plus_gammaq = 1.0 / (gammap + gammaq);
            const double gammapq = gammap * gammaq * one_o_gammap_plus_gammaq;
            const double gammap_o_gammapgammaq = gammapq * oogammaq;
            const double gammaq_o_gammapgammaq = gammapq * oogammap;
            const double Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
            const double Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
            const double Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
            const double CD_x = C[0] - D[0];
            const double CD_y = C[1] - D[1];
            const double CD_z = C[2] - D[2];
            const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

            const double PQx = Px - Qx;
            const double PQy = Py - Qy;
            const double PQz = Pz - Qz;
            const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

#if LIBINT2_DEFINED(eri, QC_x)
            const double QCx = Qx - C[0];
            erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri, QC_y)
            const double QCy = Qy - C[1];
            erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri, QC_z)
            const double QCz = Qz - C[2];
            erieval->QC_z[v] = QCz;
#endif
#if LIBINT2_DEFINED(eri, QD_x)
            const double QDx = Qx - D[0];
            erieval->QD_x[v] = QDx;
#endif
#if LIBINT2_DEFINED(eri, QD_y)
            const double QDy = Qy - D[1];
            erieval->QD_y[v] = QDy;
#endif
#if LIBINT2_DEFINED(eri, QD_z)
            const double QDz = Qz - D[2];
            erieval->QD_z[v] = QDz;
#endif

#if LIBINT2_DEFINED(eri, CD_x)
            erieval->CD_x[v] = CD_x;
#endif
#if LIBINT2_DEFINED(eri, CD_y)
            erieval->CD_y[v] = CD_y;
#endif
#if LIBINT2_DEFINED(eri, CD_z)
            erieval->CD_z[v] = CD_z;
#endif
#if LIBINT2_DEFINED(eri, DC_x)
            erieval->DC_x[v] = -CD_x;
#endif
#if LIBINT2_DEFINED(eri, DC_y)
            erieval->DC_y[v] = -CD_y;
#endif
#if LIBINT2_DEFINED(eri, DC_z)
            erieval->DC_z[v] = -CD_z;
#endif
#if LIBINT2_DEFINED(eri, oo2e)
            erieval->oo2e[v] = 0.5 * oogammaq;
#endif

            // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_x)
            erieval->TwoPRepITR_pfac0_0_0_x[v] =
                -(alpha1 * AB_x + alpha3 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_y)
            erieval->TwoPRepITR_pfac0_0_0_y[v] =
                -(alpha1 * AB_y + alpha3 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_z)
            erieval->TwoPRepITR_pfac0_0_0_z[v] =
                -(alpha1 * AB_z + alpha3 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_x)
            erieval->TwoPRepITR_pfac0_1_0_x[v] =
                -(alpha1 * AB_x + alpha3 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_y)
            erieval->TwoPRepITR_pfac0_1_0_y[v] =
                -(alpha1 * AB_y + alpha3 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_z)
            erieval->TwoPRepITR_pfac0_1_0_z[v] =
                -(alpha1 * AB_z + alpha3 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
            erieval->TwoPRepITR_pfac0_0_1_x[v] =
                (alpha0 * AB_x + alpha2 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_y)
            erieval->TwoPRepITR_pfac0_0_1_y[v] =
                (alpha0 * AB_y + alpha2 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_z)
            erieval->TwoPRepITR_pfac0_0_1_z[v] =
                (alpha0 * AB_z + alpha2 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_x)
            erieval->TwoPRepITR_pfac0_1_1_x[v] =
                (alpha0 * AB_x + alpha2 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_y)
            erieval->TwoPRepITR_pfac0_1_1_y[v] =
                (alpha0 * AB_y + alpha2 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_z)
            erieval->TwoPRepITR_pfac0_1_1_z[v] =
                (alpha0 * AB_z + alpha2 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, eoz)
            erieval->eoz[v] = gammaq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, zoe)
            erieval->zoe[v] = gammap * oogammaq;
#endif

            const double Wx =
                (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
            const double Wy =
                (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
            const double Wz =
                (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri, WP_x)
            erieval->WP_x[v] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri, WP_y)
            erieval->WP_y[v] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri, WP_z)
            erieval->WP_z[v] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri, WQ_x)
            erieval->WQ_x[v] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
            erieval->WQ_y[v] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
            erieval->WQ_z[v] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
            erieval->oo2ze[v] = 0.5 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, roz)
            erieval->roz[v] = gammapq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, roe)
            erieval->roe[v] = gammapq * oogammaq;
#endif

            if (deriv_order > 0) {
              // prefactors for derivative ERI relations
#if LIBINT2_DEFINED(eri, alpha1_rho_over_zeta2)
              erieval->alpha1_rho_over_zeta2[v] =
                  alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha2_rho_over_zeta2)
              erieval->alpha2_rho_over_zeta2[v] =
                  alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha3_rho_over_eta2)
              erieval->alpha3_rho_over_eta2[v] =
                  alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_rho_over_eta2)
              erieval->alpha4_rho_over_eta2[v] =
                  alpha3 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha1_over_zetapluseta)
              erieval->alpha1_over_zetapluseta[v] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha2_over_zetapluseta)
              erieval->alpha2_over_zetapluseta[v] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha3_over_zetapluseta)
              erieval->alpha3_over_zetapluseta[v] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_over_zetapluseta)
              erieval->alpha4_over_zetapluseta[v] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha1)
              erieval->rho12_over_alpha1[v] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
              erieval->rho12_over_alpha2[v] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha3)
              erieval->rho34_over_alpha3[v] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha4)
              erieval->rho34_over_alpha4[v] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_bra)
              erieval->two_alpha0_bra[v] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_ket)
              erieval->two_alpha0_ket[v] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_bra)
              erieval->two_alpha1_bra[v] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_ket)
              erieval->two_alpha1_ket[v] = 2.0 * alpha3;
#endif
            }

            const double K1 = exp(-rhop * AB2) * oogammap;
            const double K2 = exp(-rhoq * CD2) * oogammaq;
            const double two_times_M_PI_to_25 = 34.986836655249725693;
            double pfac =
                two_times_M_PI_to_25 * K1 * K2 * sqrt(one_o_gammap_plus_gammaq);
            pfac *= c0 * c1 * c2 * c3;

            // if veclen=1, ignore the F scratch, use the erieval directly
            if (veclen == 1 && std::is_same<double, LIBINT2_REALTYPE>::value) {
              {  // this is only used for double realtype, so nothing nefarious
                 // here, just a workaround the type system
                double* ssss_ptr =
                    reinterpret_cast<double*>(erieval->LIBINT_T_SS_EREP_SS(0));
                fmeval_chebyshev.eval(ssss_ptr, PQ2 * gammapq, amtot);
              }
              LIBINT2_REALTYPE* ssss_ptr = erieval->LIBINT_T_SS_EREP_SS(0);
              for (unsigned int l = 0; l <= amtot; ++l, ++ssss_ptr)
                *ssss_ptr = *ssss_ptr * pfac;
            } else {
              fmeval_chebyshev.eval(F, PQ2 * gammapq, amtot);
              // fmeval_taylor.eval(F,PQ2*gammapq,amtot);
              {
                LIBINT2_REALTYPE* ssss_ptr =
                    erieval->LIBINT_T_SS_EREP_SS(0) + v;
                for (unsigned int l = 0; l <= amtot; ++l, ssss_ptr += veclen)
                  *ssss_ptr = pfac * F[l];
              }
            }

          }  // end of v loop
        }
      }
    }
  }  // end of primitive loops
}

template <typename LibintEval>
void prep_libint2(LibintEval* erievals, const RandomShellSet<3u>& rsqset,
                  int norm_flag, int deriv_order = 0) {
  const uint veclen = rsqset.exp[0].size();
  const uint contrdepth = rsqset.exp[0][0].size();

  const unsigned int dummy_center =
      (LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD) ? 1 : 0;
  const double* A = &(rsqset.R[0][0]);
  const double* B = &(rsqset.R[0][0]);
  const double* C = &(rsqset.R[1][0]);
  const double* D = &(rsqset.R[2][0]);

  const uint* am = rsqset.l;
  const unsigned int amtot = am[0] + am[1] + am[2] + deriv_order;
  static double F[LIBINT_MAX_AM * 3 + 6];

  uint p012 = 0;
  for (uint p0 = 0; p0 < contrdepth; p0++) {
    for (uint p1 = 0; p1 < contrdepth; p1++) {
      for (uint p2 = 0; p2 < contrdepth; p2++, p012++) {
        LibintEval* erieval = &erievals[p012];
        erieval->veclen = veclen;
#if LIBINT2_FLOP_COUNT
        erieval->nflops = erievals[0].nflops;
#endif

        for (uint v = 0; v < veclen; v++) {
          const double alpha0 =
              (dummy_center == 0) ? 0.0 : rsqset.exp[0][v][p0];
          const double alpha1 =
              (dummy_center == 1) ? 0.0 : rsqset.exp[0][v][p0];
          const double alpha2 = rsqset.exp[1][v][p1];
          const double alpha3 = rsqset.exp[2][v][p2];

          const double c0 = (dummy_center == 0) ? 1.0 : rsqset.coef[0][v][p0];
          const double c1 = (dummy_center == 1) ? 1.0 : rsqset.coef[0][v][p0];
          const double c2 = rsqset.coef[1][v][p1];
          const double c3 = rsqset.coef[2][v][p2];

          const double gammap = alpha0 + alpha1;
          const double oogammap = 1.0 / gammap;
          const double rhop = alpha0 * alpha1 * oogammap;
          const double Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
          const double Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
          const double Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
          const double PAx = Px - A[0];
          const double PAy = Py - A[1];
          const double PAz = Pz - A[2];
          const double PBx = Px - B[0];
          const double PBy = Py - B[1];
          const double PBz = Pz - B[2];
          const double AB_x = A[0] - B[0];
          const double AB_y = A[1] - B[1];
          const double AB_z = A[2] - B[2];
          const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

#if LIBINT2_DEFINED(eri, PA_x)
          erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri, PA_y)
          erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri, PA_z)
          erieval->PA_z[v] = PAz;
#endif
#if LIBINT2_DEFINED(eri, PB_x)
          erieval->PB_x[v] = PBx;
#endif
#if LIBINT2_DEFINED(eri, PB_y)
          erieval->PB_y[v] = PBy;
#endif
#if LIBINT2_DEFINED(eri, PB_z)
          erieval->PB_z[v] = PBz;
#endif

#if LIBINT2_DEFINED(eri, AB_x)
          erieval->AB_x[v] = AB_x;
#endif
#if LIBINT2_DEFINED(eri, AB_y)
          erieval->AB_y[v] = AB_y;
#endif
#if LIBINT2_DEFINED(eri, AB_z)
          erieval->AB_z[v] = AB_z;
#endif
#if LIBINT2_DEFINED(eri, BA_x)
          erieval->BA_x[v] = -AB_x;
#endif
#if LIBINT2_DEFINED(eri, BA_y)
          erieval->BA_y[v] = -AB_y;
#endif
#if LIBINT2_DEFINED(eri, BA_z)
          erieval->BA_z[v] = -AB_z;
#endif
#if LIBINT2_DEFINED(eri, oo2z)
          erieval->oo2z[v] = 0.5 * oogammap;
#endif

          const double gammaq = alpha2 + alpha3;
          const double oogammaq = 1.0 / gammaq;
          const double rhoq = alpha2 * alpha3 * oogammaq;
          const double gammapq = gammap * gammaq / (gammap + gammaq);
          const double gammap_o_gammapgammaq = gammapq * oogammaq;
          const double gammaq_o_gammapgammaq = gammapq * oogammap;
          const double Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
          const double Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
          const double Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
          const double CD_x = C[0] - D[0];
          const double CD_y = C[1] - D[1];
          const double CD_z = C[2] - D[2];
          const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri, QC_x)
          const double QCx = Qx - C[0];
          erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri, QC_y)
          const double QCy = Qy - C[1];
          erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri, QC_z)
          const double QCz = Qz - C[2];
          erieval->QC_z[v] = QCz;
#endif
#if LIBINT2_DEFINED(eri, QD_x)
          const double QDx = Qx - D[0];
          erieval->QD_x[v] = QDx;
#endif
#if LIBINT2_DEFINED(eri, QD_y)
          const double QDy = Qy - D[1];
          erieval->QD_y[v] = QDy;
#endif
#if LIBINT2_DEFINED(eri, QD_z)
          const double QDz = Qz - D[2];
          erieval->QD_z[v] = QDz;
#endif

#if LIBINT2_DEFINED(eri, CD_x)
          erieval->CD_x[v] = CD_x;
#endif
#if LIBINT2_DEFINED(eri, CD_y)
          erieval->CD_y[v] = CD_y;
#endif
#if LIBINT2_DEFINED(eri, CD_z)
          erieval->CD_z[v] = CD_z;
#endif
#if LIBINT2_DEFINED(eri, DC_x)
          erieval->DC_x[v] = -CD_x;
#endif
#if LIBINT2_DEFINED(eri, DC_y)
          erieval->DC_y[v] = -CD_y;
#endif
#if LIBINT2_DEFINED(eri, DC_z)
          erieval->DC_z[v] = -CD_z;
#endif
#if LIBINT2_DEFINED(eri, oo2e)
          erieval->oo2e[v] = 0.5 * oogammaq;
#endif

          // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_x)
          erieval->TwoPRepITR_pfac0_0_0_x[v] =
              -(alpha1 * AB_x + alpha3 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_y)
          erieval->TwoPRepITR_pfac0_0_0_y[v] =
              -(alpha1 * AB_y + alpha3 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_z)
          erieval->TwoPRepITR_pfac0_0_0_z[v] =
              -(alpha1 * AB_z + alpha3 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_x)
          erieval->TwoPRepITR_pfac0_1_0_x[v] =
              -(alpha1 * AB_x + alpha3 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_y)
          erieval->TwoPRepITR_pfac0_1_0_y[v] =
              -(alpha1 * AB_y + alpha3 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_z)
          erieval->TwoPRepITR_pfac0_1_0_z[v] =
              -(alpha1 * AB_z + alpha3 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
          erieval->TwoPRepITR_pfac0_0_1_x[v] =
              (alpha0 * AB_x + alpha2 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_y)
          erieval->TwoPRepITR_pfac0_0_1_y[v] =
              (alpha0 * AB_y + alpha2 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_z)
          erieval->TwoPRepITR_pfac0_0_1_z[v] =
              (alpha0 * AB_z + alpha2 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_x)
          erieval->TwoPRepITR_pfac0_1_1_x[v] =
              (alpha0 * AB_x + alpha2 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_y)
          erieval->TwoPRepITR_pfac0_1_1_y[v] =
              (alpha0 * AB_y + alpha2 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_z)
          erieval->TwoPRepITR_pfac0_1_1_z[v] =
              (alpha0 * AB_z + alpha2 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, eoz)
          erieval->eoz[v] = gammaq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, zoe)
          erieval->zoe[v] = gammap * oogammaq;
#endif

          const double PQx = Px - Qx;
          const double PQy = Py - Qy;
          const double PQz = Pz - Qz;
          const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
          const double Wx =
              (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
          const double Wy =
              (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
          const double Wz =
              (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri, WP_x)
          erieval->WP_x[v] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri, WP_y)
          erieval->WP_y[v] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri, WP_z)
          erieval->WP_z[v] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri, WQ_x)
          erieval->WQ_x[v] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
          erieval->WQ_y[v] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
          erieval->WQ_z[v] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
          erieval->oo2ze[v] = 0.5 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, roz)
          erieval->roz[v] = gammapq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, roe)
          erieval->roe[v] = gammapq * oogammaq;
#endif

          // prefactors for derivative ERI relations
          if (deriv_order > 0) {
#if LIBINT2_DEFINED(eri, alpha1_rho_over_zeta2)
            erieval->alpha1_rho_over_zeta2[v] =
                alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha2_rho_over_zeta2)
            erieval->alpha2_rho_over_zeta2[v] =
                alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha3_rho_over_eta2)
            erieval->alpha3_rho_over_eta2[v] =
                alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_rho_over_eta2)
            erieval->alpha4_rho_over_eta2[v] =
                alpha3 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha1_over_zetapluseta)
            erieval->alpha1_over_zetapluseta[v] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha2_over_zetapluseta)
            erieval->alpha2_over_zetapluseta[v] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha3_over_zetapluseta)
            erieval->alpha3_over_zetapluseta[v] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_over_zetapluseta)
            erieval->alpha4_over_zetapluseta[v] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha1)
            erieval->rho12_over_alpha1[v] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
            erieval->rho12_over_alpha2[v] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha3)
            erieval->rho34_over_alpha3[v] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha4)
            erieval->rho34_over_alpha4[v] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_bra)
            erieval->two_alpha0_bra[v] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_ket)
            erieval->two_alpha0_ket[v] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_bra)
            erieval->two_alpha1_bra[v] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_ket)
            erieval->two_alpha1_ket[v] = 2.0 * alpha3;
#endif
          }

          const double K1 = exp(-rhop * AB2);
          const double K2 = exp(-rhoq * CD2);
          const double two_times_M_PI_to_25 = 34.986836655249725693;
          double pfac = two_times_M_PI_to_25 * K1 * K2 /
                        (gammap * gammaq * sqrt(gammap + gammaq));
          pfac *= c0 * c1 * c2 * c3;

          libint2::FmEval_Reference2<double>::eval(F, PQ2 * gammapq, amtot);
          // fmeval_chebyshev.eval(F,PQ2*gammapq,amtot);
          // fmeval_taylor.eval(F,PQ2*gammapq,amtot);

          // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(0))
          erieval->LIBINT_T_SS_EREP_SS(0)[v] = pfac * F[0];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(1))
          erieval->LIBINT_T_SS_EREP_SS(1)[v] = pfac * F[1];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(2))
          erieval->LIBINT_T_SS_EREP_SS(2)[v] = pfac * F[2];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(3))
          erieval->LIBINT_T_SS_EREP_SS(3)[v] = pfac * F[3];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(4))
          erieval->LIBINT_T_SS_EREP_SS(4)[v] = pfac * F[4];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(5))
          erieval->LIBINT_T_SS_EREP_SS(5)[v] = pfac * F[5];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(6))
          erieval->LIBINT_T_SS_EREP_SS(6)[v] = pfac * F[6];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(7))
          erieval->LIBINT_T_SS_EREP_SS(7)[v] = pfac * F[7];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(8))
          erieval->LIBINT_T_SS_EREP_SS(8)[v] = pfac * F[8];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(9))
          erieval->LIBINT_T_SS_EREP_SS(9)[v] = pfac * F[9];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(10))
          erieval->LIBINT_T_SS_EREP_SS(10)[v] = pfac * F[10];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(11))
          erieval->LIBINT_T_SS_EREP_SS(11)[v] = pfac * F[11];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(12))
          erieval->LIBINT_T_SS_EREP_SS(12)[v] = pfac * F[12];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(13))
          erieval->LIBINT_T_SS_EREP_SS(13)[v] = pfac * F[13];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(14))
          erieval->LIBINT_T_SS_EREP_SS(14)[v] = pfac * F[14];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(15))
          erieval->LIBINT_T_SS_EREP_SS(15)[v] = pfac * F[15];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(16))
          erieval->LIBINT_T_SS_EREP_SS(16)[v] = pfac * F[16];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(17))
          erieval->LIBINT_T_SS_EREP_SS(17)[v] = pfac * F[17];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(18))
          erieval->LIBINT_T_SS_EREP_SS(18)[v] = pfac * F[18];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(19))
          erieval->LIBINT_T_SS_EREP_SS(19)[v] = pfac * F[19];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(20))
          erieval->LIBINT_T_SS_EREP_SS(20)[v] = pfac * F[20];
#endif
        }  // end of v loop
      }
    }
  }  // end of primitive loops
}

template <typename LibintEval>
void prep_libint2(LibintEval* erievals, const RandomShellSet<2u>& rsqset,
                  int norm_flag, int deriv_order = 0) {
  const uint veclen = rsqset.exp[0].size();
  const uint contrdepth = rsqset.exp[0][0].size();

  const unsigned int dummy_center1 =
      (LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD) ? 1 : 0;
  const unsigned int dummy_center2 =
      (LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD) ? 3 : 2;
  const double* A = &(rsqset.R[0][0]);
  const double* B = &(rsqset.R[0][0]);
  const double* C = &(rsqset.R[1][0]);
  const double* D = &(rsqset.R[1][0]);

  const uint* am = rsqset.l;
  const unsigned int amtot = am[0] + am[1] + deriv_order;
  static double F[LIBINT_MAX_AM * 4 + 6];

  uint p01 = 0;
  for (uint p0 = 0; p0 < contrdepth; p0++) {
    for (uint p1 = 0; p1 < contrdepth; p1++, p01++) {
      LibintEval* erieval = &erievals[p01];
      erieval->veclen = veclen;
#if LIBINT2_FLOP_COUNT
      erieval->nflops = erievals[0].nflops;
#endif

      for (uint v = 0; v < veclen; v++) {
        const double alpha0 = (dummy_center1 == 0) ? 0.0 : rsqset.exp[0][v][p0];
        const double alpha1 = (dummy_center1 == 1) ? 0.0 : rsqset.exp[0][v][p0];
        const double alpha2 = (dummy_center2 == 2) ? 0.0 : rsqset.exp[1][v][p1];
        const double alpha3 = (dummy_center2 == 3) ? 0.0 : rsqset.exp[1][v][p1];

        const double c0 = (dummy_center1 == 0) ? 1.0 : rsqset.coef[0][v][p0];
        const double c1 = (dummy_center1 == 1) ? 1.0 : rsqset.coef[0][v][p0];
        const double c2 = (dummy_center2 == 2) ? 1.0 : rsqset.coef[1][v][p1];
        const double c3 = (dummy_center2 == 3) ? 1.0 : rsqset.coef[1][v][p1];

        const double gammap = alpha0 + alpha1;
        const double oogammap = 1.0 / gammap;
        const double rhop = alpha0 * alpha1 * oogammap;
        const double Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
        const double Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
        const double Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
        const double PAx = Px - A[0];
        const double PAy = Py - A[1];
        const double PAz = Pz - A[2];
        const double PBx = Px - B[0];
        const double PBy = Py - B[1];
        const double PBz = Pz - B[2];
        const double AB_x = A[0] - B[0];
        const double AB_y = A[1] - B[1];
        const double AB_z = A[2] - B[2];
        const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

#if LIBINT2_DEFINED(eri, PA_x)
        erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri, PA_y)
        erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri, PA_z)
        erieval->PA_z[v] = PAz;
#endif
#if LIBINT2_DEFINED(eri, PB_x)
        erieval->PB_x[v] = PBx;
#endif
#if LIBINT2_DEFINED(eri, PB_y)
        erieval->PB_y[v] = PBy;
#endif
#if LIBINT2_DEFINED(eri, PB_z)
        erieval->PB_z[v] = PBz;
#endif

#if LIBINT2_DEFINED(eri, AB_x)
        erieval->AB_x[v] = AB_x;
#endif
#if LIBINT2_DEFINED(eri, AB_y)
        erieval->AB_y[v] = AB_y;
#endif
#if LIBINT2_DEFINED(eri, AB_z)
        erieval->AB_z[v] = AB_z;
#endif
#if LIBINT2_DEFINED(eri, oo2z)
        erieval->oo2z[v] = 0.5 * oogammap;
#endif

        const double gammaq = alpha2 + alpha3;
        const double oogammaq = 1.0 / gammaq;
        const double rhoq = alpha2 * alpha3 * oogammaq;
        const double gammapq = gammap * gammaq / (gammap + gammaq);
        const double gammap_o_gammapgammaq = gammapq * oogammaq;
        const double gammaq_o_gammapgammaq = gammapq * oogammap;
        const double Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
        const double Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
        const double Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
        const double CD_x = C[0] - D[0];
        const double CD_y = C[1] - D[1];
        const double CD_z = C[2] - D[2];
        const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri, QC_x)
        const double QCx = Qx - C[0];
        erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri, QC_y)
        const double QCy = Qy - C[1];
        erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri, QC_z)
        const double QCz = Qz - C[2];
        erieval->QC_z[v] = QCz;
#endif
#if LIBINT2_DEFINED(eri, QD_x)
        const double QDx = Qx - D[0];
        erieval->QD_x[v] = QDx;
#endif
#if LIBINT2_DEFINED(eri, QD_y)
        const double QDy = Qy - D[1];
        erieval->QD_y[v] = QDy;
#endif
#if LIBINT2_DEFINED(eri, QD_z)
        const double QDz = Qz - D[2];
        erieval->QD_z[v] = QDz;
#endif

#if LIBINT2_DEFINED(eri, CD_x)
        erieval->CD_x[v] = CD_x;
#endif
#if LIBINT2_DEFINED(eri, CD_y)
        erieval->CD_y[v] = CD_y;
#endif
#if LIBINT2_DEFINED(eri, CD_z)
        erieval->CD_z[v] = CD_z;
#endif
#if LIBINT2_DEFINED(eri, oo2e)
        erieval->oo2e[v] = 0.5 * oogammaq;
#endif

        // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_x)
        erieval->TwoPRepITR_pfac0_0_0_x[v] =
            -(alpha1 * AB_x + alpha3 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_y)
        erieval->TwoPRepITR_pfac0_0_0_y[v] =
            -(alpha1 * AB_y + alpha3 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_z)
        erieval->TwoPRepITR_pfac0_0_0_z[v] =
            -(alpha1 * AB_z + alpha3 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_x)
        erieval->TwoPRepITR_pfac0_1_0_x[v] =
            -(alpha1 * AB_x + alpha3 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_y)
        erieval->TwoPRepITR_pfac0_1_0_y[v] =
            -(alpha1 * AB_y + alpha3 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_z)
        erieval->TwoPRepITR_pfac0_1_0_z[v] =
            -(alpha1 * AB_z + alpha3 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
        erieval->TwoPRepITR_pfac0_0_1_x[v] =
            (alpha0 * AB_x + alpha2 * CD_x) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_y)
        erieval->TwoPRepITR_pfac0_0_1_y[v] =
            (alpha0 * AB_y + alpha2 * CD_y) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_z)
        erieval->TwoPRepITR_pfac0_0_1_z[v] =
            (alpha0 * AB_z + alpha2 * CD_z) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_x)
        erieval->TwoPRepITR_pfac0_1_1_x[v] =
            (alpha0 * AB_x + alpha2 * CD_x) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_y)
        erieval->TwoPRepITR_pfac0_1_1_y[v] =
            (alpha0 * AB_y + alpha2 * CD_y) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_z)
        erieval->TwoPRepITR_pfac0_1_1_z[v] =
            (alpha0 * AB_z + alpha2 * CD_z) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, eoz)
        erieval->eoz[v] = gammaq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, zoe)
        erieval->zoe[v] = gammap * oogammaq;
#endif

        const double PQx = Px - Qx;
        const double PQy = Py - Qy;
        const double PQz = Pz - Qz;
        const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
        const double Wx =
            (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
        const double Wy =
            (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
        const double Wz =
            (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri, WP_x)
        erieval->WP_x[v] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri, WP_y)
        erieval->WP_y[v] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri, WP_z)
        erieval->WP_z[v] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri, WQ_x)
        erieval->WQ_x[v] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
        erieval->WQ_y[v] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
        erieval->WQ_z[v] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
        erieval->oo2ze[v] = 0.5 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, roz)
        erieval->roz[v] = gammapq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, roe)
        erieval->roe[v] = gammapq * oogammaq;
#endif

        // prefactors for derivative ERI relations
        if (deriv_order > 0) {
#if LIBINT2_DEFINED(eri, alpha1_rho_over_zeta2)
          erieval->alpha1_rho_over_zeta2[v] =
              alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha2_rho_over_zeta2)
          erieval->alpha2_rho_over_zeta2[v] =
              alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri, alpha3_rho_over_eta2)
          erieval->alpha3_rho_over_eta2[v] =
              alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_rho_over_eta2)
          erieval->alpha4_rho_over_eta2[v] =
              alpha3 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha1_over_zetapluseta)
          erieval->alpha1_over_zetapluseta[v] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha2_over_zetapluseta)
          erieval->alpha2_over_zetapluseta[v] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha3_over_zetapluseta)
          erieval->alpha3_over_zetapluseta[v] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_over_zetapluseta)
          erieval->alpha4_over_zetapluseta[v] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha1)
          erieval->rho12_over_alpha1[v] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
          erieval->rho12_over_alpha2[v] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha3)
          erieval->rho34_over_alpha3[v] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha4)
          erieval->rho34_over_alpha4[v] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_bra)
          erieval->two_alpha0_bra[v] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_ket)
          erieval->two_alpha0_ket[v] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_bra)
          erieval->two_alpha1_bra[v] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_ket)
          erieval->two_alpha1_ket[v] = 2.0 * alpha3;
#endif
        }

        const double K1 = exp(-rhop * AB2);
        const double K2 = exp(-rhoq * CD2);
        const double two_times_M_PI_to_25 = 34.986836655249725693;
        double pfac = two_times_M_PI_to_25 * K1 * K2 /
                      (gammap * gammaq * sqrt(gammap + gammaq));
        pfac *= c0 * c1 * c2 * c3;

        libint2::FmEval_Reference2<double>::eval(F, PQ2 * gammapq, amtot);
        // fmeval_chebyshev.eval(F,PQ2*gammapq,amtot);
        // fmeval_taylor.eval(F,PQ2*gammapq,amtot);

        // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(0))
        erieval->LIBINT_T_SS_EREP_SS(0)[v] = pfac * F[0];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(1))
        erieval->LIBINT_T_SS_EREP_SS(1)[v] = pfac * F[1];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(2))
        erieval->LIBINT_T_SS_EREP_SS(2)[v] = pfac * F[2];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(3))
        erieval->LIBINT_T_SS_EREP_SS(3)[v] = pfac * F[3];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(4))
        erieval->LIBINT_T_SS_EREP_SS(4)[v] = pfac * F[4];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(5))
        erieval->LIBINT_T_SS_EREP_SS(5)[v] = pfac * F[5];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(6))
        erieval->LIBINT_T_SS_EREP_SS(6)[v] = pfac * F[6];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(7))
        erieval->LIBINT_T_SS_EREP_SS(7)[v] = pfac * F[7];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(8))
        erieval->LIBINT_T_SS_EREP_SS(8)[v] = pfac * F[8];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(9))
        erieval->LIBINT_T_SS_EREP_SS(9)[v] = pfac * F[9];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(10))
        erieval->LIBINT_T_SS_EREP_SS(10)[v] = pfac * F[10];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(11))
        erieval->LIBINT_T_SS_EREP_SS(11)[v] = pfac * F[11];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(12))
        erieval->LIBINT_T_SS_EREP_SS(12)[v] = pfac * F[12];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(13))
        erieval->LIBINT_T_SS_EREP_SS(13)[v] = pfac * F[13];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(14))
        erieval->LIBINT_T_SS_EREP_SS(14)[v] = pfac * F[14];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(15))
        erieval->LIBINT_T_SS_EREP_SS(15)[v] = pfac * F[15];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(16))
        erieval->LIBINT_T_SS_EREP_SS(16)[v] = pfac * F[16];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(17))
        erieval->LIBINT_T_SS_EREP_SS(17)[v] = pfac * F[17];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(18))
        erieval->LIBINT_T_SS_EREP_SS(18)[v] = pfac * F[18];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(19))
        erieval->LIBINT_T_SS_EREP_SS(19)[v] = pfac * F[19];
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(20))
        erieval->LIBINT_T_SS_EREP_SS(20)[v] = pfac * F[20];
#endif
      }  // end of v loop
    }
  }  // end of primitive loops
}
