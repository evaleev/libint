//////////////////////////////////////////////////////////////////////////////////
// This program shows how to use Libint for computing electron repulsion integrals
// (c) 2011 Edward F. Valeev
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <libint2.h>

// angular momentum labels
const char cglabel[] = { 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n',
                         'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' };

#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

/// this function computes all data required by Libint to evaluate the integrals efficiently
template<typename LibintEval>
void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3], int norm_flag);

/*!
 eri()

 This is a very inefficient function for computing ERIs
 over primitive Gaussian functions. The argument
 list is self-explanatory, except for norm_flag:

 \param norm_flag:  tells what kind of normalization to use,
 0 - no normalization, >0 - normalized ERI
 */

double eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1,
           const double* A, unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, const double* B, unsigned int l3, unsigned int m3,
           unsigned int n3, double alpha3, const double* C, unsigned int l4,
           unsigned int m4, unsigned int n4, double alpha4, const double* D,
           int norm_flag);

int main(int argc, char** argv) {

  typedef unsigned int uint;

  // this initializes internal Libint data structures -- must happen once in the program
  LIBINT2_PREFIXED_NAME(libint2_static_init)();

  // This example assumes that your library does not support for vectorization
  const uint veclen = 1;

  // These parameters define 4 primitive Gaussian basis functions
  double alpha[4] = { 0.5, 1.0, 1.5, 2.0 }; // orbital exponents for the 4 functions
  double A[3] = { 1.0, 2.0, 3.0 }; // position of function 1
  double B[3] = { 1.5, 2.5, 3.5 }; // etc.
  double C[3] = { 4.0, 2.0, 0.0 };
  double D[3] = { 3.0, 3.0, 1.0 };

  // LIBINT2_MAX_AM_ERI is a macro defined in libint2.h that specifies the maximum angular momentum
  // this Libint library instance can handle
  // const unsigned int ammax = LIBINT2_MAX_AM_ERI;
  const unsigned int ammax = 3;
  // Libint_t is the type of a data structure used to pass basis function data to Libint
  Libint_t inteval;
  // Libint_t objects must be initialized prior to use
  // your program may have several such objects to implement computation of integrals in multiple threads
  // each thread would have its own Libint_t object
  LIBINT2_PREFIXED_NAME( libint2_init_eri)(&inteval, ammax, 0);

  // I will evaluate every type of integrals supported by a given library
  // change lmax above to evaluate a subset of integrals
  for (uint am0 = 0; am0 <= ammax; ++am0) {
    for (uint am1 = 0; am1 <= ammax; ++am1) {
      for (uint am2 = 0; am2 <= ammax; ++am2) {
        for (uint am3 = 0; am3 <= ammax; ++am3) {

          // Libint only provides code to computes a symmetry-unique subset of integrals
          // the default convention is to compute integrals (0 1|2 3) with am0 >= am1, am2 >= am3, and am2+am3 >= am0+am1
          // in general you will need to permute shells to satisfy this ordering
          // most production codes do this by precomputing "significant" shell pair lists
          bool can_compute = (am0 >= am1) && (am2 >= am3) &&
                             (am0 + am1 <= am2 + am3);
          can_compute &= (am0 != 0 || am1 != 0 || am2 != 0 || am3 != 0); // skip (ss|ss) integral -- no need to use Libint for that one
          if (can_compute == false)
            continue;

          // this function fills in all data expected by Libint to evaluate the given integral set
          prep_libint2(&inteval, am0, alpha[0], A, am1, alpha[1], B, am2,
                       alpha[2], C, am3, alpha[3], D, 0);

          // announce what we are computing now
          std::cout << "Testing (" << cglabel[am0] << cglabel[am1] << "|"
              << cglabel[am2] << cglabel[am3] << ") ... ";
          std::cout.flush();

          // these are uncontracted basis functions
          // N.B. to compute integrals over contracted functions allocate an array of Libint_t objects (one for each primitive combination)
          //      and fill each with primitive combination data. You only need to call libint2_init_eri using pointer to the first object!
          inteval.contrdepth = 1;
          // this compute the shell set (quartet) of integrals
          LIBINT2_PREFIXED_NAME
              ( libint2_build_eri)[am0][am1][am2][am3](&inteval);

          bool success = true;
          int ijkl = 0;
          // verify each integral against a painfully slow but bulletproof ERI function
          // iterate over basis functions in shell 0
          // this doubly-nested loop iterates over functions in Libint standard ordering; for example, for f functions is:
          // xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
          // Libint compiler can produce libraries that support different orderings
          for (uint k0 = 0; k0 <= am0; ++k0) {
            for (uint n0 = 0; n0 <= k0; ++n0) {
              const uint m0 = k0 - n0;
              const uint l0 = am0 - k0;
              // l0, m0, n0 are x,y,z exponents for the first basis function

              for (uint k1 = 0; k1 <= am1; ++k1) {
                for (uint n1 = 0; n1 <= k1; ++n1) {
                  const uint m1 = k1 - n1;
                  const uint l1 = am1 - k1;

                  for (uint k2 = 0; k2 <= am2; ++k2) {
                    for (uint n2 = 0; n2 <= k2; ++n2) {
                      const uint m2 = k2 - n2;
                      const uint l2 = am2 - k2;

                      for (uint k3 = 0; k3 <= am3; ++k3) {
                        for (uint n3 = 0; n3 <= k3; ++n3) {
                          const uint m3 = k3 - n3;
                          const uint l3 = am3 - k3;

                          // eri() computes reference value for the ERI over 4 primitive functions
                          // it is extremely slow and is only used here for validation of Libint output
                          const double ref_eri =
                              eri(l0, m0, n0, alpha[0], A, l1, m1, n1,
                                  alpha[1], B, l2, m2, n2, alpha[2], C, l3, m3,
                                  n3, alpha[3], D, 0);

                          //
                          const double libint_eri = inteval.targets[0][ijkl];

                          if (fabs((ref_eri - libint_eri) / libint_eri)
                              > 1.0E-8) {
#if 1
                            std::cout << std::endl << "Elem " << ijkl
                            << " : eri.cc = " << ref_eri << " libint = "
                            << libint_eri;
#endif
                            success = false;
                          }

                          ++ijkl;

                        }
                      }
                    }
                  }
                }
              }
            }
          } // end of loop over basis functions in the shell quartet

          std::cout << (success ? "ok" : "failed") << std::endl;

        }
      }
    }
  } // end of loop over angular momenta

  // this releases all memory that was allocated for this object
  LIBINT2_PREFIXED_NAME( libint2_cleanup_eri)(&inteval);

  return 0;
}

/////////
// eri() computes electron repulsion integrals using very slow but reliable method
/////////

static double *df;
static double *fac;
static double **bc;

void calc_f(double *, int, double);
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                  double alpha1, const double* A);
double* init_array(unsigned long int size);
double** block_matrix(unsigned long int nrow, unsigned long int ncol);
void free_array(double* array);

double eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1,
           const double* A, unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, const double* B, unsigned int l3, unsigned int m3,
           unsigned int n3, double alpha3, const double* C, unsigned int l4,
           unsigned int m4, unsigned int n4, double alpha4, const double* D,
           int norm_flag) {

  const double gammap = alpha1 + alpha2;
  const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
  const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
  const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
  const double PAx = Px - A[0];
  const double PAy = Py - A[1];
  const double PAz = Pz - A[2];
  const double PBx = Px - B[0];
  const double PBy = Py - B[1];
  const double PBz = Pz - B[2];
  const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1]
      - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);

  const double gammaq = alpha3 + alpha4;
  const double gammapq = gammap * gammaq / (gammap + gammaq);
  const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
  const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
  const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
  const double QCx = Qx - C[0];
  const double QCy = Qy - C[1];
  const double QCz = Qz - C[2];
  const double QDx = Qx - D[0];
  const double QDy = Qy - D[1];
  const double QDz = Qz - D[2];
  const double CD2 = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1]
      - D[1]) + (C[2] - D[2]) * (C[2] - D[2]);

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

  int u1, u2, v1, v2, w1, w2, tx, ty, tz, txmax, tymax, tzmax;
  int i, j, k;
  int lp, lq, mp, mq, np, nq;
  int zeta;
  double *flp, *flq, *fmp, *fmq, *fnp, *fnq;
  double *F;
  double K1, K2;
  double Gx, Gy, Gz;
  double pfac;
  double result = 0.0;
  double tmp;
  int u1max, u2max, v1max, v2max, w1max, w2max;

  K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
  K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
  pfac = 2 * std::pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
      * sqrt(gammap + gammaq));

  if (fac == NULL) {
    fac = init_array(MAXFAC);
    fac[0] = 1.0;
    for (i = 1; i < MAXFAC; i++)
      fac[i] = fac[i - 1] * i;
    bc = block_matrix(MAXFAC, MAXFAC);
    for (i = 0; i < MAXFAC; i++)
      for (j = 0; j <= i; j++)
        bc[i][j] = fac[i] / (fac[i - j] * fac[j]);
  }

  if (norm_flag > 0) {
    pfac *= norm_const(l1, m1, n1, alpha1, A);
    pfac *= norm_const(l2, m2, n2, alpha2, B);
    pfac *= norm_const(l3, m3, n3, alpha3, C);
    pfac *= norm_const(l4, m4, n4, alpha4, D);
  }

  F = init_array(l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4 + 1);
  calc_f(F, l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4,
         PQ2 * gammapq);

  flp = init_array(l1 + l2 + 1);
  for (k = 0; k <= l1 + l2; k++)
    for (i = 0; i <= MIN(k,l1); i++) {
      j = k - i;
      if (j > l2)
        continue;
      tmp = bc[l1][i] * bc[l2][j];
      if (l1 - i > 0)
        tmp *= pow(PAx, l1 - i);
      if (l2 - j > 0)
        tmp *= pow(PBx, l2 - j);
      flp[k] += tmp;
    }
  fmp = init_array(m1 + m2 + 1);
  for (k = 0; k <= m1 + m2; k++)
    for (i = 0; i <= MIN(k,m1); i++) {
      j = k - i;
      if (j > m2)
        continue;
      tmp = bc[m1][i] * bc[m2][j];
      if (m1 - i > 0)
        tmp *= pow(PAy, m1 - i);
      if (m2 - j > 0)
        tmp *= pow(PBy, m2 - j);
      fmp[k] += tmp;
    }
  fnp = init_array(n1 + n2 + 1);
  for (k = 0; k <= n1 + n2; k++)
    for (i = 0; i <= MIN(k,n1); i++) {
      j = k - i;
      if (j > n2)
        continue;
      tmp = bc[n1][i] * bc[n2][j];
      if (n1 - i > 0)
        tmp *= pow(PAz, n1 - i);
      if (n2 - j > 0)
        tmp *= pow(PBz, n2 - j);
      fnp[k] += tmp;
    }
  flq = init_array(l3 + l4 + 1);
  for (k = 0; k <= l3 + l4; k++)
    for (i = 0; i <= MIN(k,l3); i++) {
      j = k - i;
      if (j > l4)
        continue;
      tmp = bc[l3][i] * bc[l4][j];
      if (l3 - i > 0)
        tmp *= pow(QCx, l3 - i);
      if (l4 - j > 0)
        tmp *= pow(QDx, l4 - j);
      flq[k] += tmp;
    }
  fmq = init_array(m3 + m4 + 1);
  for (k = 0; k <= m3 + m4; k++)
    for (i = 0; i <= MIN(k,m3); i++) {
      j = k - i;
      if (j > m4)
        continue;
      tmp = bc[m3][i] * bc[m4][j];
      if (m3 - i > 0)
        tmp *= pow(QCy, m3 - i);
      if (m4 - j > 0)
        tmp *= pow(QDy, m4 - j);
      fmq[k] += tmp;
    }
  fnq = init_array(n3 + n4 + 1);
  for (k = 0; k <= n3 + n4; k++)
    for (i = 0; i <= MIN(k,n3); i++) {
      j = k - i;
      if (j > n4)
        continue;
      tmp = bc[n3][i] * bc[n4][j];
      if (n3 - i > 0)
        tmp *= pow(QCz, n3 - i);
      if (n4 - j > 0)
        tmp *= pow(QDz, n4 - j);
      fnq[k] += tmp;
    }

  for (lp = 0; lp <= l1 + l2; lp++)
    for (lq = 0; lq <= l3 + l4; lq++) {
      u1max = lp / 2;
      u2max = lq / 2;
      for (u1 = 0; u1 <= u1max; u1++)
        for (u2 = 0; u2 <= u2max; u2++) {
          Gx = pow(-1, lp) * flp[lp] * flq[lq] * fac[lp] * fac[lq]
              * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq) * fac[lp + lq - 2
              * u1 - 2 * u2] * pow(gammapq, lp + lq - 2 * u1 - 2 * u2)
              / (fac[u1] * fac[u2] * fac[lp - 2 * u1] * fac[lq - 2 * u2]);
          for (mp = 0; mp <= m1 + m2; mp++)
            for (mq = 0; mq <= m3 + m4; mq++) {
              v1max = mp / 2;
              v2max = mq / 2;
              for (v1 = 0; v1 <= v1max; v1++)
                for (v2 = 0; v2 <= v2max; v2++) {
                  Gy = pow(-1, mp) * fmp[mp] * fmq[mq] * fac[mp] * fac[mq]
                      * pow(gammap, v1 - mp) * pow(gammaq, v2 - mq) * fac[mp
                      + mq - 2 * v1 - 2 * v2] * pow(gammapq,
                                                    mp + mq - 2 * v1 - 2 * v2)
                      / (fac[v1] * fac[v2] * fac[mp - 2 * v1]
                          * fac[mq - 2 * v2]);
                  for (np = 0; np <= n1 + n2; np++)
                    for (nq = 0; nq <= n3 + n4; nq++) {
                      w1max = np / 2;
                      w2max = nq / 2;
                      for (w1 = 0; w1 <= w1max; w1++)
                        for (w2 = 0; w2 <= w2max; w2++) {
                          Gz = pow(-1, np) * fnp[np] * fnq[nq] * fac[np]
                              * fac[nq] * pow(gammap, w1 - np) * pow(gammaq,
                                                                     w2 - nq)
                              * fac[np + nq - 2 * w1 - 2 * w2]
                              * pow(gammapq, np + nq - 2 * w1 - 2 * w2)
                              / (fac[w1] * fac[w2] * fac[np - 2 * w1] * fac[nq
                                  - 2 * w2]);
                          txmax = (lp + lq - 2 * u1 - 2 * u2) / 2;
                          tymax = (mp + mq - 2 * v1 - 2 * v2) / 2;
                          tzmax = (np + nq - 2 * w1 - 2 * w2) / 2;
                          for (tx = 0; tx <= txmax; tx++)
                            for (ty = 0; ty <= tymax; ty++)
                              for (tz = 0; tz <= tzmax; tz++) {
                                zeta = lp + lq + mp + mq + np + nq - 2 * u1 - 2
                                    * u2 - 2 * v1 - 2 * v2 - 2 * w1 - 2 * w2
                                    - tx - ty - tz;
                                result += Gx * Gy * Gz * F[zeta]
                                    * pow(-1, tx + ty + tz) * pow(
                                                                  PQx,
                                                                  lp + lq - 2
                                                                      * u1 - 2
                                                                      * u2 - 2
                                                                      * tx)
                                    * pow(PQy,
                                          mp + mq - 2 * v1 - 2 * v2 - 2 * ty)
                                    * pow(PQz,
                                          np + nq - 2 * w1 - 2 * w2 - 2 * tz)
                                    / (pow(
                                           4,
                                           u1 + u2 + tx + v1 + v2 + ty + w1
                                               + w2 + tz) * pow(gammapq, tx)
                                        * pow(gammapq, ty) * pow(gammapq, tz)
                                        * fac[lp + lq - 2 * u1 - 2 * u2 - 2
                                            * tx] * fac[tx] * fac[mp + mq - 2
                                        * v1 - 2 * v2 - 2 * ty] * fac[ty]
                                        * fac[np + nq - 2 * w1 - 2 * w2 - 2
                                            * tz] * fac[tz]);
                              }
                        }
                    }
                }
            }
        }
    }

  free_array(F);
  free_array(flp);
  free_array(fmp);
  free_array(fnp);
  free_array(flq);
  free_array(fmq);
  free_array(fnq);

  return result * pfac;
}

/*!
 calc_f()

 This function computes infamous integral Fn(t). For its definition
 see Obara and Saika paper, or Shavitt's chapter in the
 Methods in Computational Physics book (see reference below).
 This piece of code is from Dr. Justin Fermann's program CINTS

 \ingroup (QT)
 */
void calc_f(double *F, int n, double t) {
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0 / M_2_SQRTPI;
  double et;

  if (df == NULL) {
    df = init_array(2 * MAXFAC);
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (i = 3; i < MAXFAC * 2; i++) {
      df[i] = (i - 1) * df[i - 2];
    }
  }

  if (t > 20.0) { /* For big t's do upward recursion */
    t2 = 2 * t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K * erf(t) / t;
    for (m = 0; m <= n - 1; m++) {
      F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
    }
  } else { /* For smaller t's compute F with highest n using
   asymptotic series (see I. Shavitt in
   Methods in Computational Physics, ed. B. Alder eta l,
   vol 2, 1963, page 8) */
    et = exp(-t);
    t2 = 2 * t;
    m2 = 2 * n;
    num = df[m2];
    i = 0;
    sum = 1.0 / (m2 + 1);
    do {
      i++;
      num = num * t2;
      term1 = num / df[m2 + 2 * i + 2];
      sum += term1;
    } while (fabs(term1) > EPS && i < MAXFAC);
    F[n] = sum * et;
    for (m = n - 1; m >= 0; m--) { /* And then do downward recursion */
      F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
    }
  }
}

/*!
 norm_const()

 \ingroup (QT)
 */
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                  double alpha1, const double* A) {
  return pow(2 * alpha1 / M_PI, 0.75) * pow(4 * alpha1, 0.5 * (l1 + m1 + n1))
      / sqrt(df[2 * l1] * df[2 * m1] * df[2 * n1]);
}

double* init_array(unsigned long int size) {
  double* result = new double[size];
  for (int i = 0; i < size; i++)
    result[i] = 0.0;
  return result;
}

double** block_matrix(unsigned long int nrow, unsigned long int ncol) {
  double** rows = new double*[nrow];
  rows[0] = new double[nrow * ncol];
  for (int i = 1; i < nrow; i++)
    rows[i] = rows[i - 1] + ncol;

  return rows;
}

void free_array(double* array) {
  delete[] array;
}

template<typename LibintEval>
void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3], int norm_flag) {

  const unsigned int am = am1 + am2 + am3 + am4;
  double* F = init_array(am + 1);

  const double gammap = alpha1 + alpha2;
  const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
  const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
  const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
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

#if LIBINT2_DEFINED(eri,PA_x)
  erieval->PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
  erieval->PA_y[0] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
  erieval->PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
  erieval->PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
  erieval->PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
  erieval->PB_z[0] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
  erieval->AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
  erieval->AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
  erieval->AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
  erieval->oo2z[0] = 0.5/gammap;
#endif

  const double gammaq = alpha3 + alpha4;
  const double gammapq = gammap * gammaq / (gammap + gammaq);
  const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
  const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
  const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
  const double QCx = Qx - C[0];
  const double QCy = Qy - C[1];
  const double QCz = Qz - C[2];
  const double QDx = Qx - D[0];
  const double QDy = Qy - D[1];
  const double QDz = Qz - D[2];
  const double CD_x = C[0] - D[0];
  const double CD_y = C[1] - D[1];
  const double CD_z = C[2] - D[2];
  const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri,QC_x)
  erieval->QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
  erieval->QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
  erieval->QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
  erieval->QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
  erieval->QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
  erieval->QD_z[0] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
  erieval->CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
  erieval->CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
  erieval->CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
  erieval->oo2e[0] = 0.5/gammaq;
#endif

  // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
  erieval->TwoPRepITR_pfac0_0_x[0] = - (alpha1*AB_x + alpha3*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
  erieval->TwoPRepITR_pfac0_0_y[0] = - (alpha1*AB_y + alpha3*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
  erieval->TwoPRepITR_pfac0_0_z[0] = - (alpha1*AB_z + alpha3*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
  erieval->TwoPRepITR_pfac0_1_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
  erieval->TwoPRepITR_pfac0_1_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
  erieval->TwoPRepITR_pfac0_1_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_0)
  erieval->TwoPRepITR_pfac1_0[0] = -gammaq/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_1)
  erieval->TwoPRepITR_pfac1_1[0] = -gammap/gammaq;
#endif

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
  const double Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
  const double Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
  const double Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

#if LIBINT2_DEFINED(eri,WP_x)
  erieval->WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
  erieval->WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
  erieval->WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
  erieval->WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
  erieval->WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
  erieval->WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
  erieval->oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
  erieval->roz[0] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
  erieval->roe[0] = gammapq/gammaq;
#endif

  double K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
  double K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
  double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
      * sqrt(gammap + gammaq));

  if (norm_flag > 0) {
    /*    pfac *= norm_const(l1,m1,n1,alpha1,A);
     pfac *= norm_const(l2,m2,n2,alpha2,B);
     pfac *= norm_const(l3,m3,n3,alpha3,C);
     pfac *= norm_const(l4,m4,n4,alpha4,D);*/
  }

  calc_f(F, am, PQ2 * gammapq);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20];
#endif

  free_array(F);
}
