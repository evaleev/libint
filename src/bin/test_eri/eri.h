/*
 *  Copyright (C) 2004-2017 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_testeri_eri_h_
#define _libint2_src_bin_testeri_eri_h_

#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdlib.h>

#include <libint2/boys.h>

#if defined(__cplusplus)
#if HAVE_MPFR
# include <cstddef>
# include <gmpxx.h>
# include <mpfr.h>
 typedef mpf_class LIBINT2_REF_REALTYPE;
 /// implement exp for mpf_class using MPFR ... I do not claim to know what issues the rounding presents here
 inline mpf_class exp(mpf_class x) {
   mpfr_t x_r; mpfr_init(x_r);
   mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

   mpfr_t expx_r; mpfr_init(expx_r);
   mpfr_exp(expx_r, x_r, MPFR_RNDN);

   mpf_t expx;
   mpf_init(expx);
   mpfr_get_f(expx, expx_r, MPFR_RNDN);
   mpf_class result(expx);
   return result;
 }
 /// implement pow for mpf_class using MPFR ... I do not claim to know what issues the rounding presents here
 inline mpf_class pow(mpf_class x, int a) {
   mpf_t x_to_a;
   mpf_init(x_to_a);
   if (a >= 0)
     mpf_pow_ui(x_to_a, x.get_mpf_t(), (unsigned int) a);
   else
     mpf_pow_ui(x_to_a, x.get_mpf_t(), (unsigned int) (-a));
   mpf_class result(x_to_a);
   if (a < 0)
     result = 1.0 / result;
   return result;
 }
 /// implement erf for mpf_class using MPFR ... I do not claim to know what issues the rounding presents here
 inline mpf_class erf(mpf_class x) {
   mpfr_t x_r; mpfr_init(x_r);
   mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

   mpfr_t erfx_r; mpfr_init(erfx_r);
   mpfr_erf(erfx_r, x_r, MPFR_RNDN);

   mpf_t erfx;
   mpf_init(erfx);
   mpfr_get_f(erfx, erfx_r, MPFR_RNDN);
   mpf_class result(erfx);
   return result;
 }
 /// implement log for mpf_class using MPFR ... I do not claim to know what issues the rounding presents here
 inline mpf_class log(mpf_class x) {
   mpfr_t x_r; mpfr_init(x_r);
   mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

   mpfr_t logx_r; mpfr_init(logx_r);
   mpfr_log(logx_r, x_r, MPFR_RNDN);

   mpf_t logx;
   mpf_init(logx);
   mpfr_get_f(logx, logx_r, MPFR_RNDN);
   mpf_class result(logx);
   return result;
 }

#else
 typedef double LIBINT2_REF_REALTYPE;
#endif
#else
 typedef double LIBINT2_REF_REALTYPE;
#endif

#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */

struct ExpensiveMath {
#if defined(__clusplus)
    typedef mpz_class BigInt;
#else
    typedef double BigInt;
#endif
    BigInt *df;
    BigInt *fac;
    BigInt **bc;
    const static unsigned int MAXFAC = 100;

    ExpensiveMath() {
      df = new BigInt[2*MAXFAC];
      df[0] = 1ul;
      df[1] = 1ul;
      df[2] = 1ul;
      for(unsigned long i=3; i<MAXFAC*2; i++) {
        df[i] = (i-1)*df[i-2];
      }

      fac = new BigInt[MAXFAC];
      fac[0] = 1ul;
      for(unsigned long i=1;i<MAXFAC;i++)
        fac[i] = fac[i-1]*i;

      bc = new BigInt*[MAXFAC];
      bc[0] = new BigInt[MAXFAC*MAXFAC];
      for(unsigned long i=1;i<MAXFAC;i++)
        bc[i] = bc[i-1] + MAXFAC;
      for(unsigned long i=0;i<MAXFAC;i++) {
        bc[i][0] = 1ul;
        for(unsigned long j=1;j<=i;j++) {
          bc[i][j] = (bc[i][j-1] * (i-j+1) )/ j; // fac[i]/(fac[i-j]*fac[j]);
        }
      }
    }

    ~ExpensiveMath() {
      delete[] df;
      delete[] fac;
      delete[] bc[0];
      delete[] bc;
    }

    template <typename Real>
    Real norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                    Real alpha1, const Real* A)
    {
      const Real sqrt_twoalpha1_over_pi = sqrt(2*alpha1/M_PI);
      const Real sqrtsqrt_twoalpha1_over_pi = sqrt(sqrt_twoalpha1_over_pi);
      const Real sqrt_alpha1 = sqrt(alpha1);

      return sqrt_twoalpha1_over_pi * sqrtsqrt_twoalpha1_over_pi * (2 * pow(sqrt_alpha1,l1+m1+n1))/sqrt(df[2*l1]*df[2*m1]*df[2*n1]);
    }
};

/*!
  calc_f()

  This function computes infamous integral Fn(t). For its definition
  see Obara and Saika paper, or Shavitt's chapter in the
  Methods in Computational Physics book (see reference below).
  This piece of code is from Dr. Justin Fermann's program CINTS

 \ingroup (QT)
*/
template <typename Real> void calc_f(Real* F, int n, Real t)
{
  static ExpensiveMath expmath;
  int i, m, k;
  int m2;
  Real t2;
  Real num;
  Real sum;
  Real term1, term2;
  static Real K = 1.0/M_2_SQRTPI;
  Real et;

  if (t>20.0){   /* For big t's do upward recursion */
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
    }
  }
  else {        /* For smaller t's compute F with highest n using
                   asymptotic series (see I. Shavitt in
                   Methods in Computational Physics, ed. B. Alder eta l,
                   vol 2, 1963, page 8) */
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = expmath.df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/expmath.df[m2+2*i+2];
      sum += term1;
    } while (term1 > EPS && i < expmath.MAXFAC);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){        /* And then do downward recursion */
      F[m] = (t2*F[m+1] + et)/(2*m+1);
    }
  }
}

namespace {
  template <typename Int>
  int parity(Int a) {
    if (a%2 == 1) return -1;
    return 1;
  }
}

namespace libint2 {
  int min(int t1, unsigned int t2) { return t1 > (int)t2 ? (int)t2 : t1; }
};

/*!
  eri()

  This is a very inefficient function for computing ERIs
  over primitive Gaussian functions. The argument
  list is self-explanatory, except for norm_flag:

  \param norm_flag:  tells what kind of normalization to use,
         0 - no normalization, >0 - normalized ERI
*/
inline LIBINT2_REF_REALTYPE eri(unsigned int l1, unsigned int m1, unsigned int n1,
                                LIBINT2_REF_REALTYPE alpha1, const LIBINT2_REF_REALTYPE* A,
                                unsigned int l2, unsigned int m2, unsigned int n2,
                                LIBINT2_REF_REALTYPE alpha2, const LIBINT2_REF_REALTYPE* B,
                                unsigned int l3, unsigned int m3, unsigned int n3,
                                LIBINT2_REF_REALTYPE alpha3, const LIBINT2_REF_REALTYPE* C,
                                unsigned int l4, unsigned int m4, unsigned int n4,
                                LIBINT2_REF_REALTYPE alpha4, const LIBINT2_REF_REALTYPE* D,
                                int norm_flag)
{
  static ExpensiveMath expmath;

  const LIBINT2_REF_REALTYPE gammap = alpha1 + alpha2;
  const LIBINT2_REF_REALTYPE Px = (alpha1*A[0] + alpha2*B[0])/gammap;
  const LIBINT2_REF_REALTYPE Py = (alpha1*A[1] + alpha2*B[1])/gammap;
  const LIBINT2_REF_REALTYPE Pz = (alpha1*A[2] + alpha2*B[2])/gammap;
  const LIBINT2_REF_REALTYPE PAx = Px - A[0];
  const LIBINT2_REF_REALTYPE PAy = Py - A[1];
  const LIBINT2_REF_REALTYPE PAz = Pz - A[2];
  const LIBINT2_REF_REALTYPE PBx = Px - B[0];
  const LIBINT2_REF_REALTYPE PBy = Py - B[1];
  const LIBINT2_REF_REALTYPE PBz = Pz - B[2];
  const LIBINT2_REF_REALTYPE AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);

  const LIBINT2_REF_REALTYPE gammaq = alpha3 + alpha4;
  const LIBINT2_REF_REALTYPE gammapq = gammap*gammaq/(gammap+gammaq);
  const LIBINT2_REF_REALTYPE Qx = (alpha3*C[0] + alpha4*D[0])/gammaq;
  const LIBINT2_REF_REALTYPE Qy = (alpha3*C[1] + alpha4*D[1])/gammaq;
  const LIBINT2_REF_REALTYPE Qz = (alpha3*C[2] + alpha4*D[2])/gammaq;
  const LIBINT2_REF_REALTYPE QCx = Qx - C[0];
  const LIBINT2_REF_REALTYPE QCy = Qy - C[1];
  const LIBINT2_REF_REALTYPE QCz = Qz - C[2];
  const LIBINT2_REF_REALTYPE QDx = Qx - D[0];
  const LIBINT2_REF_REALTYPE QDy = Qy - D[1];
  const LIBINT2_REF_REALTYPE QDz = Qz - D[2];
  const LIBINT2_REF_REALTYPE CD2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);

  const LIBINT2_REF_REALTYPE PQx = Px - Qx;
  const LIBINT2_REF_REALTYPE PQy = Py - Qy;
  const LIBINT2_REF_REALTYPE PQz = Pz - Qz;
  const LIBINT2_REF_REALTYPE PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;

  const LIBINT2_REF_REALTYPE four = 4;

  int u1,u2,v1,v2,w1,w2,tx,ty,tz,txmax,tymax,tzmax;
  int i,j,k;
  int lp,lq,mp,mq,np,nq;
  int zeta;
  LIBINT2_REF_REALTYPE *flp, *flq, *fmp, *fmq, *fnp, *fnq;
  LIBINT2_REF_REALTYPE *F;
  LIBINT2_REF_REALTYPE K1, K2;
  LIBINT2_REF_REALTYPE Gx,Gy,Gz;
  LIBINT2_REF_REALTYPE pfac;
  LIBINT2_REF_REALTYPE result = 0.0;
  LIBINT2_REF_REALTYPE tmp;
  int u1max,u2max,v1max,v2max,w1max,w2max;

  K1 = exp(-alpha1*alpha2*AB2/gammap);
  K2 = exp(-alpha3*alpha4*CD2/gammaq);
  const LIBINT2_REF_REALTYPE m_pi = M_PI;
  const LIBINT2_REF_REALTYPE m_pi_2 = m_pi * m_pi;
  pfac = 2*sqrt(m_pi_2 * m_pi_2 * m_pi)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));

  if (norm_flag > 0) {
    pfac *= expmath.norm_const(l1,m1,n1,alpha1,A);
    pfac *= expmath.norm_const(l2,m2,n2,alpha2,B);
    pfac *= expmath.norm_const(l3,m3,n3,alpha3,C);
    pfac *= expmath.norm_const(l4,m4,n4,alpha4,D);
  }

  const int ltot = l1+l2+l3+l4+m1+m2+m3+m4+n1+n2+n3+n4;
  F = new LIBINT2_REF_REALTYPE[ltot+1];
  //calc_f<LIBINT2_REF_REALTYPE>(F,ltot,PQ2*gammapq);
  libint2::FmEval_Reference2<LIBINT2_REF_REALTYPE>::eval(F,PQ2*gammapq,ltot,1e-15);

#define DEBUG 0
#if DEBUG
  std::cout << "pfac = " << pfac << "  K1 = " << K1 << "  K2 = " << K2 << std::endl;
  std::cout << "PQ2 = " << PQ2 << "  gammapq = " << gammapq << std::endl;
  for(int i=0; i<=ltot; ++i)
    std::cout << "F[" << i << "] = " << F[i] << std::endl;
#endif

  flp = new LIBINT2_REF_REALTYPE[l1+l2+1];
  for(k=0;k<=l1+l2;k++) {
    flp[k] = 0.0;
    for(i=0;i<=libint2::min(k,l1);i++) {
      j = k-i;
      if (j > l2) continue;
      tmp = expmath.bc[l1][i]*expmath.bc[l2][j];
      if (l1 - i > 0)
        tmp *= pow(PAx, l1 - i);
      if (l2 - j > 0)
        tmp *= pow(PBx, l2 - j);
      flp[k] += tmp;
    }
  }
  fmp = new LIBINT2_REF_REALTYPE[m1+m2+1];
  for(k=0;k<=m1+m2;k++) {
    fmp[k] = 0.0;
    for(i=0;i<=libint2::min(k,m1);i++) {
      j = k-i;
      if (j > m2) continue;
      tmp = expmath.bc[m1][i]*expmath.bc[m2][j];
      if (m1 - i > 0)
        tmp *= pow(PAy, m1 - i);
      if (m2 - j > 0)
        tmp *= pow(PBy, m2 - j);
      fmp[k] += tmp;
    }
  }
  fnp = new LIBINT2_REF_REALTYPE[n1+n2+1];
  for(k=0;k<=n1+n2;k++) {
    fnp[k] = 0.0;
    for(i=0;i<=libint2::min(k,n1);i++) {
      j = k-i;
      if (j > n2) continue;
      tmp = expmath.bc[n1][i]*expmath.bc[n2][j];
      if (n1 - i > 0)
        tmp *= pow(PAz, n1 - i);
      if (n2 - j > 0)
        tmp *= pow(PBz, n2 - j);
      fnp[k] += tmp;
    }
  }
  flq = new LIBINT2_REF_REALTYPE[l3+l4+1];
  for(k=0;k<=l3+l4;k++) {
    flq[k] = 0.0;
    for(i=0;i<=libint2::min(k,l3);i++) {
      j = k-i;
      if (j > l4) continue;
      tmp = expmath.bc[l3][i]*expmath.bc[l4][j];
      if (l3 - i > 0)
        tmp *= pow(QCx, l3 - i);
      if (l4 - j > 0)
        tmp *= pow(QDx, l4 - j);
      flq[k] += tmp;
    }
  }
  fmq = new LIBINT2_REF_REALTYPE[m3+m4+1];
  for(k=0;k<=m3+m4;k++) {
    fmq[k] = 0.0;
    for(i=0;i<=libint2::min(k,m3);i++) {
      j = k-i;
      if (j > m4) continue;
      tmp = expmath.bc[m3][i]*expmath.bc[m4][j];
      if (m3 - i > 0)
        tmp *= pow(QCy, m3 - i);
      if (m4 - j > 0)
        tmp *= pow(QDy, m4 - j);
      fmq[k] += tmp;
    }
  }
  fnq = new LIBINT2_REF_REALTYPE[n3+n4+1];
  for(k=0;k<=n3+n4;k++) {
    fnq[k] = 0.0;
    for(i=0;i<=libint2::    min(k,n3);i++) {
      j = k-i;
      if (j > n4) continue;
      tmp = expmath.bc[n3][i]*expmath.bc[n4][j];
      if (n3 - i > 0)
        tmp *= pow(QCz, n3 - i);
      if (n4 - j > 0)
        tmp *= pow(QDz, n4 - j);
      fnq[k] += tmp;
    }
  }

  for (lp = 0; lp <= l1 + l2; lp++) {
    for (lq = 0; lq <= l3 + l4; lq++) {
      u1max = lp / 2;
      u2max = lq / 2;
      for (u1 = 0; u1 <= u1max; u1++) {
        for (u2 = 0; u2 <= u2max; u2++) {
          Gx = parity(lp) * flp[lp] * flq[lq] * expmath.fac[lp]
              * expmath.fac[lq] * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq)
              * expmath.fac[lp + lq - 2 * u1 - 2 * u2]
              * pow(gammapq, lp + lq - 2 * u1 - 2 * u2)
              / (expmath.fac[u1] * expmath.fac[u2] * expmath.fac[lp - 2 * u1]
                  * expmath.fac[lq - 2 * u2]);
          for (mp = 0; mp <= m1 + m2; mp++) {
            for (mq = 0; mq <= m3 + m4; mq++) {
              v1max = mp / 2;
              v2max = mq / 2;
              for (v1 = 0; v1 <= v1max; v1++) {
                for (v2 = 0; v2 <= v2max; v2++) {
                  Gy =
                      parity(mp) * fmp[mp] * fmq[mq] * expmath.fac[mp]
                          * expmath.fac[mq] * pow(gammap, v1 - mp)
                          * pow(gammaq, v2 - mq)
                          * expmath.fac[mp + mq - 2 * v1 - 2 * v2]
                          * pow(gammapq, mp + mq - 2 * v1 - 2 * v2)
                          / (expmath.fac[v1] * expmath.fac[v2]
                              * expmath.fac[mp - 2 * v1]
                              * expmath.fac[mq - 2 * v2]);
                  for (np = 0; np <= n1 + n2; np++) {
                    for (nq = 0; nq <= n3 + n4; nq++) {
                      w1max = np / 2;
                      w2max = nq / 2;
                      for (w1 = 0; w1 <= w1max; w1++) {
                        for (w2 = 0; w2 <= w2max; w2++) {
                          Gz = parity(np) * fnp[np] * fnq[nq] * expmath.fac[np]
                              * expmath.fac[nq] * pow(gammap, w1 - np)
                              * pow(gammaq, w2 - nq)
                              * expmath.fac[np + nq - 2 * w1 - 2 * w2]
                              * pow(gammapq, np + nq - 2 * w1 - 2 * w2)
                              / (expmath.fac[w1] * expmath.fac[w2]
                                  * expmath.fac[np - 2 * w1]
                                  * expmath.fac[nq - 2 * w2]);
                          txmax = (lp + lq - 2 * u1 - 2 * u2) / 2;
                          tymax = (mp + mq - 2 * v1 - 2 * v2) / 2;
                          tzmax = (np + nq - 2 * w1 - 2 * w2) / 2;
                          for (tx = 0; tx <= txmax; tx++) {
                            for (ty = 0; ty <= tymax; ty++) {
                              for (tz = 0; tz <= tzmax; tz++) {
                                zeta = lp + lq + mp + mq + np + nq - 2 * u1
                                    - 2 * u2 - 2 * v1 - 2 * v2 - 2 * w1 - 2 * w2
                                    - tx - ty - tz;
                                result += Gx * Gy * Gz * F[zeta]
                                    * parity(tx + ty + tz)
                                    * pow(PQx,
                                          lp + lq - 2 * u1 - 2 * u2 - 2 * tx)
                                    * pow(PQy,
                                          mp + mq - 2 * v1 - 2 * v2 - 2 * ty)
                                    * pow(PQz,
                                          np + nq - 2 * w1 - 2 * w2 - 2 * tz)
                                    / (pow(
                                        four,
                                        u1 + u2 + tx + v1 + v2 + ty + w1 + w2
                                            + tz) * pow(gammapq, tx)
                                        * pow(gammapq, ty) * pow(gammapq, tz)
                                        * expmath.fac[lp + lq - 2 * u1 - 2 * u2
                                            - 2 * tx] * expmath.fac[tx]
                                        * expmath.fac[mp + mq - 2 * v1 - 2 * v2
                                            - 2 * ty] * expmath.fac[ty]
                                        * expmath.fac[np + nq - 2 * w1 - 2 * w2
                                            - 2 * tz] * expmath.fac[tz]);
#if DEBUG
                                std::cout << "Gx = " << Gx << std::endl;
                                std::cout << "Gy = " << Gy << std::endl;
                                std::cout << "Gz = " << Gz << std::endl;
                                std::cout << "result = " << result << std::endl;
#endif
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] F;
  delete[] flp;
  delete[] fmp;
  delete[] fnp;
  delete[] flq;
  delete[] fmq;
  delete[] fnq;

  return result*pfac;
}


/// same as above, except specifies derivative order; uses the above function
/// to compute
/// \param deriv_order a set of 12 integers (3 coordinates for each of 4 basis function origins)
inline LIBINT2_REF_REALTYPE eri(const unsigned int* deriv_index,
                                unsigned int l1, unsigned int m1, unsigned int n1,
                                LIBINT2_REF_REALTYPE alpha1, const LIBINT2_REF_REALTYPE* A,
                                unsigned int l2, unsigned int m2, unsigned int n2,
                                LIBINT2_REF_REALTYPE alpha2, const LIBINT2_REF_REALTYPE* B,
                                unsigned int l3, unsigned int m3, unsigned int n3,
                                LIBINT2_REF_REALTYPE alpha3, const LIBINT2_REF_REALTYPE* C,
                                unsigned int l4, unsigned int m4, unsigned int n4,
                                LIBINT2_REF_REALTYPE alpha4, const LIBINT2_REF_REALTYPE* D,
                                int norm_flag) {
  const unsigned int deriv_order = std::accumulate(deriv_index, deriv_index+12, 0u);
  if (deriv_order == 0)
    return eri(l1, m1, n1, alpha1, A,
               l2, m2, n2, alpha2, B,
               l3, m3, n3, alpha3, C,
               l4, m4, n4, alpha4, D,
               norm_flag);
  else {

    struct nonzero_t {
        bool operator()(unsigned int val) {
          return val > 0;
        }
    };

    nonzero_t nonzero;
    // handle one derivative at a time
    const unsigned int* di_ptr = std::find_if(deriv_index, deriv_index+12, nonzero);
    const unsigned int di = di_ptr - deriv_index;
    const unsigned int dc = di / 3;
    const unsigned int dxyz = di % 3;

    unsigned int qn[4][3];
    qn[0][0] = l1;
    qn[0][1] = m1;
    qn[0][2] = n1;
    qn[1][0] = l2;
    qn[1][1] = m2;
    qn[1][2] = n2;
    qn[2][0] = l3;
    qn[2][1] = m3;
    qn[2][2] = n3;
    qn[3][0] = l4;
    qn[3][1] = m4;
    qn[3][2] = n4;

    LIBINT2_REF_REALTYPE alpha[4];
    alpha[0] = alpha1;
    alpha[1] = alpha2;
    alpha[2] = alpha3;
    alpha[3] = alpha4;

    ++qn[dc][dxyz];
    unsigned int deriv_index_minus1[12];  std::copy(deriv_index, deriv_index+12, deriv_index_minus1);
    --deriv_index_minus1[di];
    LIBINT2_REF_REALTYPE result = eri(deriv_index_minus1,
                        qn[0][0], qn[0][1], qn[0][2], alpha[0], A,
                        qn[1][0], qn[1][1], qn[1][2], alpha[1], B,
                        qn[2][0], qn[2][1], qn[2][2], alpha[2], C,
                        qn[3][0], qn[3][1], qn[3][2], alpha[3], D,
                        norm_flag) *
                    2.0 * alpha[dc];
    --qn[dc][dxyz];

    if (qn[dc][dxyz] > 0) {
      --qn[dc][dxyz];
      result -= eri(deriv_index_minus1,
                    qn[0][0], qn[0][1], qn[0][2], alpha[0], A,
                    qn[1][0], qn[1][1], qn[1][2], alpha[1], B,
                    qn[2][0], qn[2][1], qn[2][2], alpha[2], C,
                    qn[3][0], qn[3][1], qn[3][2], alpha[3], D,
                    norm_flag) *
                    (qn[dc][dxyz] + 1);
      ++qn[dc][dxyz];
    }
    return result;
  }
  assert(false); // unreachable;
}




#endif // header guard
