/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_bin_testeri_eri_h_
#define _libint2_src_bin_testeri_eri_h_

#include <libint2/config.h>

#if defined(__cplusplus)
#if HAVE_MPFR
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
#else
 typedef double LIBINT2_REF_REALTYPE;
#endif
#else
 typedef double LIBINT2_REF_REALTYPE;
#endif

     /*!
       eri()

       This is a very inefficient function for computing ERIs
       over primitive Gaussian functions. The argument
       list is self-explanatory, except for norm_flag:

       \param norm_flag:  tells what kind of normalization to use,
              0 - no normalization, >0 - normalized ERI
     */
LIBINT2_REF_REALTYPE eri(unsigned int l1, unsigned int m1, unsigned int n1,
                         LIBINT2_REF_REALTYPE alpha1, const LIBINT2_REF_REALTYPE* A,
                         unsigned int l2, unsigned int m2, unsigned int n2,
                         LIBINT2_REF_REALTYPE alpha2, const LIBINT2_REF_REALTYPE* B,
                         unsigned int l3, unsigned int m3, unsigned int n3,
                         LIBINT2_REF_REALTYPE alpha3, const LIBINT2_REF_REALTYPE* C,
                         unsigned int l4, unsigned int m4, unsigned int n4,
                         LIBINT2_REF_REALTYPE alpha4, const LIBINT2_REF_REALTYPE* D, int norm_flag);

/// same as above, except specifies derivative order; uses the above function
/// to compute
/// \param deriv_order a set of 12 integers (3 coordinates for each of 4 basis function origins)
LIBINT2_REF_REALTYPE eri(const unsigned int* deriv_index,
                         unsigned int l1, unsigned int m1, unsigned int n1,
                         LIBINT2_REF_REALTYPE alpha1, const LIBINT2_REF_REALTYPE* A,
                         unsigned int l2, unsigned int m2, unsigned int n2,
                         LIBINT2_REF_REALTYPE alpha2, const LIBINT2_REF_REALTYPE* B,
                         unsigned int l3, unsigned int m3, unsigned int n3,
                         LIBINT2_REF_REALTYPE alpha3, const LIBINT2_REF_REALTYPE* C,
                         unsigned int l4, unsigned int m4, unsigned int n4,
                         LIBINT2_REF_REALTYPE alpha4, const LIBINT2_REF_REALTYPE* D, int norm_flag);

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

    template <typename Real>
    Real norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                    Real alpha1, const Real* A)
    {
      return pow(2*alpha1/M_PI,0.75)*pow(4*alpha1,0.5*(l1+m1+n1))/sqrt(df[2*l1]*df[2*m1]*df[2*n1]);
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

#if 0
namespace libint2 {

  class ExpensiveMath {
    public:
      ExpensiveMath(int ifac, int idf);
      ~ExpensiveMath();
      double *fac;
      double *df;
  };

  #define TAYLOR_INTERPOLATION_AND_RECURSION 0  // compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only

  /**
   * Computes the Boys function using Taylor expansion in (T-T_k), where T_k is the closest value of T that exists in
   * the table of precomputed values.
   *
   * @tparam RealP the real number data type used to precompute Fm(T). The default is double. Must implement the following operators and functions:
   *    \li arithmetic (+, -, *, /)
   *    \li log, exp, sqrt, pow, floor
   * @tparam RealE the real number data type used to evaluate Fm(T). The default is same as RealP. Must implement the following operators and functions:
   *    \li arithmetic (+, -, *, /)
   *    \li log, exp, sqrt, pow, floor
   * @tparam InterpolationOrder the interpolation order. The default is 6.
   */
  template <typename RealP = double, typename RealE = RealP, unsigned int InterpolationOrder = 6u>
  class FmTaylor {
    public:
      static const unsigned int max_interp_order = 8u;

      /// Initialize to compute Fm with 0<=m<=\c m_max, and estimated relative precision \c prec
      FmTaylor(unsigned int m_max, Real prec);
      ~FmTaylor();
      /// returns the pointer to the result buffer
      Real *values(unsigned int J, Real T);

    private:
      static Real relative_zero_;
      Real **grid_; /* Table of "exact" Fm(T) values. Row index corresponds to
       values of T (max_T+1 rows), column index to values
       of m (max_m+1 columns) */
      Real delT_; /* The step size for T, depends on cutoff */
      Real oodelT_; /* 1.0 / delT_, see above */
      Real cutoff_; /* Tolerance cutoff used in all computations of Fj(T) */
      unsigned int max_m_; /* Maximum value of m in the table, depends on cutoff
       and the number of terms in Taylor interpolation */
      unsigned int max_T_; /* Maximum index of T in the table, depends on cutoff
       and m */
      Real *T_crit_; /* Maximum T for each row, depends on cutoff;
       for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
       for a given m and T_idx > max_T_idx[m] use the asymptotic formula */
      Real *F_; /* Here computed values of Fj(T) are stored */

      ExpensiveMath<Real> ExpMath_;
  };

}; // namespace libint2
#endif

#endif // header guard
