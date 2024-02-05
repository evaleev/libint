/*
 *  Copyright (C) 2004-2024 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <assert.h>
#include <libint2.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

Libint_t erieval;
double* F;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/* computes F[0] .. F[max_m] */
extern void calc_f(double* F, double T, unsigned int max_m);

#if defined(LIBINT2_SUPPORT_ERI) && LIBINT2_MAX_AM_eri >= 1

/** This function evaluates ERI over 4 primitive Gaussian shells.
    See tests/eri/test.cc for an example of how to deal with
    contracted Gaussians.

    For simplicity, many details are omitted here, e.g. normalization.
  */
void _compute_eri(Libint_t* erieval, unsigned int am1, double alpha1, double* A,
                  unsigned int am2, double alpha2, double* B, unsigned int am3,
                  double alpha3, double* C, unsigned int am4, double alpha4,
                  double* D) {
  /* I will assume that libint2_static_init() and
   * libint2_init_eri(&erieval,max_am,0) had been called elsewhere! */

  double gammap, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, AB2;
  double gammaq, Qx, Qy, Qz, QCx, QCy, QCz, QDx, QDy, QDz, CD2;
  double gammapq, PQx, PQy, PQz, PQ2, Wx, Wy, Wz;
  double K1, K2, pfac;
  unsigned int am;
  double* eri_shell_set;

  /*
     Compute requisite data -- many of these quantities would be precomputed
     for all nonnegligible shell pairs somewhere else
  */
  gammap = alpha1 + alpha2;
  Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
  Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
  Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
  PAx = Px - A[0];
  PAy = Py - A[1];
  PAz = Pz - A[2];
  PBx = Px - B[0];
  PBy = Py - B[1];
  PBz = Pz - B[2];
  AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) +
        (A[2] - B[2]) * (A[2] - B[2]);

  erieval->PA_x[0] = PAx;
  erieval->PA_y[0] = PAy;
  erieval->PA_z[0] = PAz;
  erieval->AB_x[0] = A[0] - B[0];
  erieval->AB_y[0] = A[1] - B[1];
  erieval->AB_z[0] = A[2] - B[2];
  erieval->oo2z[0] = 0.5 / gammap;

  gammaq = alpha3 + alpha4;
  gammapq = gammap * gammaq / (gammap + gammaq);
  Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
  Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
  Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
  QCx = Qx - C[0];
  QCy = Qy - C[1];
  QCz = Qz - C[2];
  QDx = Qx - D[0];
  QDy = Qy - D[1];
  QDz = Qz - D[2];
  CD2 = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1] - D[1]) +
        (C[2] - D[2]) * (C[2] - D[2]);

  erieval->QC_x[0] = QCx;
  erieval->QC_y[0] = QCy;
  erieval->QC_z[0] = QCz;
  erieval->CD_x[0] = C[0] - D[0];
  erieval->CD_y[0] = C[1] - D[1];
  erieval->CD_z[0] = C[2] - D[2];
  erieval->oo2e[0] = 0.5 / gammaq;

  PQx = Px - Qx;
  PQy = Py - Qy;
  PQz = Pz - Qz;
  PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
  Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
  Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
  Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

  erieval->WP_x[0] = Wx - Px;
  erieval->WP_y[0] = Wy - Py;
  erieval->WP_z[0] = Wz - Pz;
  erieval->WQ_x[0] = Wx - Qx;
  erieval->WQ_y[0] = Wy - Qy;
  erieval->WQ_z[0] = Wz - Qz;
  erieval->oo2ze[0] = 0.5 / (gammap + gammaq);
  erieval->roz[0] = gammapq / gammap;
  erieval->roe[0] = gammapq / gammaq;

  K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
  K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
  pfac =
      2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq * sqrt(gammap + gammaq));

  /*
     evaluate Boys function F_m for all m in [0,am]
  */
  am = am1 + am2 + am3 + am4;
  calc_f(F, PQ2 * gammapq, am);

  /* (00|00)^m = pfac * F_m */
  assert(am <= 4);
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac * F[0];
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac * F[1];
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac * F[2];
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac * F[3];
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac * F[4];

  /* compute ERIs */
  libint2_build_eri[am1][am2][am3][am4](erieval);
}

void init_c_api(unsigned int max_am) {
  libint2_init_eri(&erieval, max_am, 0);
  F = malloc(sizeof(double) * (4 * max_am + 1));
#if LIBINT_CONTRACTED_INTS
  /* if have support for contracted integrals, set the contraction length to 1
   */
  erieval.contrdepth = 1;
#endif
}

double* compute_eri(unsigned int am1, double alpha1, double* A,
                    unsigned int am2, double alpha2, double* B,
                    unsigned int am3, double alpha3, double* C,
                    unsigned int am4, double alpha4, double* D) {
  _compute_eri(&erieval, am1, alpha1, A, am2, alpha2, B, am3, alpha3, C, am4,
               alpha4, D);
  return erieval.targets[0];
}

void finalize_c_api() {
  free(F);
  libint2_cleanup_eri(&erieval);
}

#endif
