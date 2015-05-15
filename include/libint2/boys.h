/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

// prototype for the Boys function engines (Boys function = Fm(T))
// the Chebyshev extrapolation code is based on that by Frank Neese

#ifndef _libint2_src_lib_libint_boys_h_
#define _libint2_src_lib_libint_boys_h_

#if defined(__cplusplus)

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <libint2/vector.h>
#include <cassert>
#include <vector>
#include <algorithm>

// some features require at least C++11
#if __cplusplus > 199711L
#include <memory>
#endif

#if HAVE_LAPACK // use F77-type interface for now, switch to LAPACKE later
extern "C" void dgesv_(const int* n,
                       const int* nrhs, double* A, const int* lda,
                       int* ipiv, double* b, const int* ldb,
                       int* info);
#endif

namespace libint2 {

  /// holds tables of expensive quantities
  template<typename Real>
  class ExpensiveNumbers {
    public:
      ExpensiveNumbers(int ifac = -1, int idf = -1, int ibc = -1) {
        if (ifac >= 0) {
          fac.resize(ifac + 1);
          fac[0] = 1.0;
          for (int i = 1; i <= ifac; i++) {
            fac[i] = i * fac[i - 1];
          }
        }

        if (idf >= 0) {
          df.resize(idf + 1);
          /* df[i] gives (i-1)!!, so that (-1)!! is defined... */
          df[0] = 1.0;
          if (idf >= 1)
            df[1] = 1.0;
          if (idf >= 2)
            df[2] = 1.0;
          for (int i = 3; i <= idf; i++) {
            df[i] = (i - 1) * df[i - 2];
          }
        }

        if (ibc >= 0) {
          bc_.resize((ibc+1)*(ibc+1));
          std::fill(bc_.begin(), bc_.end(), Real(0));
          bc.resize(ibc+1);
          bc[0] = &bc_[0];
          for(int i=1; i<=ibc; ++i)
            bc[i] = bc[i-1] + (ibc+1);

          for(int i=0; i<=ibc; i++)
            bc[i][0] = 1.0;
          for(int i=0; i<=ibc; i++)
            for(int j=1; j<=i; ++j)
              bc[i][j] = bc[i][j-1] * Real(i-j+1) / Real(j);
        }

        for (int i = 0; i < 128; i++) {
          twoi1[i] = 1.0 / (Real(2.0) * i + Real(1.0));
          ihalf[i] = Real(i) - Real(0.5);
        }

      }

      ~ExpensiveNumbers() {
      }

      std::vector<Real> fac;
      std::vector<Real> df;
      std::vector<Real*> bc;

      // these quantitites are needed with indices <= mmax
      // 64 is sufficient to handle up to 4 center integrals with up to L=15 basis functions
      // but need higher values for Yukawa integrals ...
      Real twoi1[128]; /* 1/(2 i + 1); needed for downward recursion */
      Real ihalf[128]; /* i - 0.5, needed for upward recursion */

    private:
      std::vector<Real> bc_;
  };

#define _local_min_macro(a,b) ((a) > (b) ? (a) : (b))

  /** Computes the Boys function, \f$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \, {\rm d}u \f$,
    * using single algorithm (asymptotic expansion). Slow for the sake of precision control.
    * Useful in two cases:
    * <ul>
    *   <li> for reference purposes, if \c Real supports high/arbitrary precision, and </li>
    *   <li> for moderate values of \f$ T \f$, if \c Real is a low-precision floating-point type.
    *        N.B. FmEval_Reference2 , which can compute for all practical values of \f$ T \f$ and \f$ m \f$, is recommended
    *        with standard \c Real types (\c double and \c float). </li>
    * </ul>
    *
    * \note Precision is controlled heuristically, i.e. cannot be guaranteed mathematically;
    *       will stop if absolute precision is reached, or precision of \c Real is exhausted.
    *       It is important that \c std::numeric_limits<Real> is defined appropriately.
    *
    * @tparam Real the type to use for all floating-point computations.
    *         Must be able to compute logarithm and exponential, i.e.
    *         log(x) and exp(x), where x is Real, must be valid expressions.
    */
  template<typename Real>
  struct FmEval_Reference {

      /// computes a single value of \f$ F_m(T) \f$ using MacLaurin series.
      static Real eval(Real T, size_t m, Real absolute_precision) {
        assert(m < 100);
        static const Real T_crit = std::numeric_limits<Real>::is_bounded == true ? -log( std::numeric_limits<Real>::min() * 100.5 / 2. ) : Real(0) ;
        if (std::numeric_limits<Real>::is_bounded && T > T_crit)
          throw std::overflow_error("FmEval_Reference<Real>::eval: Real lacks precision for the given value of argument T");
        Real denom = (m + 0.5);
        Real term = 0.5 * exp(-T) / denom;
        Real old_term = 0.0;
        Real sum = term;
        //Real rel_error;
        Real epsilon;
        const Real relative_zero = std::numeric_limits<Real>::epsilon();
        const Real absolute_precision_o_1000 = absolute_precision * 0.001;
        do {
          denom += 1.0;
          old_term = term;
          term = old_term * T / denom;
          sum += term;
          //rel_error = term / sum;
          // stop if adding a term smaller or equal to absolute_precision/1000 and smaller than relative_zero * sum
          // When sum is small in absolute value, the second threshold is more important
          epsilon = _local_min_macro(absolute_precision_o_1000, sum*relative_zero);
        } while (term > epsilon || old_term < term);

        return sum;
      }

      /// fills up an array of Fm(T) for m in [0,mmax]
      /// @param[out] Fm array to be filled in with the Boys function values, must be at least mmax+1 elements long
      /// @param[in] T the Boys function argument
      /// @param[in] mmax the maximum value of m for which Boys function will be computed;
      /// @param[in] absolute_precision the absolute precision to which to compute the result
      static void eval(Real* Fm, Real T, size_t mmax, Real absolute_precision) {

        // evaluate for mmax using MacLaurin series
        // it converges fastest for the largest m -> use it to compute Fmmax(T)
        //  see JPC 94, 5564 (1990).
        for(size_t m=0; m<=mmax; ++m)
          Fm[m] = eval(T, m, absolute_precision);
        return;
        /** downward recursion does not maintain absolute precision, only relative precision, and cannot be used for T > 10
        if (mmax > 0) {
          const Real T2 = 2.0 * T;
          const Real exp_T = exp(-T);
          for (int m = mmax - 1; m >= 0; m--)
            Fm[m] = (Fm[m + 1] * T2 + exp_T) / (2 * m + 1);
        }
        */
      }

  };

  /** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \, {\rm d}u \$,
    * using multi-algorithm approach (upward precision for T>=30, and asymptotic summation for T<30).
    * This is slow and should be used for reference purposes, e.g. computing the interpolation tables.
    * Precision is not always guaranteed as it is limited by the precision of \c Real type.
    * When \c Real is \c double, can maintain 1e-14 precision for up to m=38 and 0<=T<=1e9 .
    *
    * @tparam Real the type to use for all floating-point computations.
    *         Must be able to compute logarithm, exponential, square root, and error function, i.e.
    *         log(x), exp(x), sqrt(x), and erf(x), where x is Real, must be valid expressions.
    */
  template<typename Real>
  struct FmEval_Reference2 {

      /// fills up an array of Fm(T) for m in [0,mmax]
      /// @param[out] Fm array to be filled in with the Boys function values, must be at least mmax+1 elements long
      /// @param[in] t the Boys function argument
      /// @param[in] mmax the maximum value of m for which Boys function will be computed;
      /// @param[in] absolute_precision the absolute precision to which to compute the result
      static void eval(Real* Fm, Real t, size_t mmax, Real absolute_precision) {

        if (t < Real(30)) {
          FmEval_Reference<Real>::eval(Fm,t,mmax,absolute_precision);
        }
        else {
          const Real two_over_sqrt_pi{1.128379167095512573896158903121545171688101258657997713688171443421284936882986828973487320404214727};
          const Real K = 1.0/two_over_sqrt_pi;

          auto t2 = 2*t;
          auto et = exp(-t);
          auto sqrt_t = sqrt(t);
          Fm[0] = K*erf(sqrt_t)/sqrt_t;
          if (mmax > 0)
          for(size_t m=0; m<=mmax-1; m++) {
            Fm[m+1] = ((2*m + 1)*Fm[m] - et)/(t2);
          }
        }
      }

  };

  /** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \, {\rm d}u \$,
    * using Chebyshev interpolation.
    * based on the code from ORCA by Dr. Frank Neese.
    */
  template <typename Real = double>
  class FmEval_Chebyshev3 {

      static const int NGRID = 4096; //!< number of grid points
      static const int INTERPOLATION_ORDER = 4;   //!< interpolation order + 1
      static const bool INTERPOLATION_AND_RECURSION = false; //!< compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only

      const Real T_crit;          //!< critical value of T above which safe to use upward recusion
      const Real delta;           //!< grid size
      const Real one_over_delta;  //! 1/delta
      int mmax;                     //!< the maximum m that is tabulated
      ExpensiveNumbers<double> numbers_;
      Real *c; /* the Chebyshev coefficients table, NGRID by mmax*interpolation_order */

    public:
      /// \param m_max maximum value of the Boys function index; set to -1 to skip initialization
      /// \param precision the desired precision
      FmEval_Chebyshev3(int m_max, double = 0.0) :
          T_crit(30.0), // this translates in appr. 1e-15  error in upward recursion, see the note below
          delta(T_crit / (NGRID - 1)),
          one_over_delta(1.0 / delta),
          mmax(m_max), numbers_(14, 0) {
        assert(mmax <= 63);
        if (m_max >= 0)
          init();
      }
      ~FmEval_Chebyshev3() {
        if (mmax >= 0) {
          free(c);
        }
      }

      // some features require at least C++11
#if __cplusplus > 199711L
      /// Singleton interface allows to manage the lone instance; adjusts max m values as needed in thread-safe fashion
      static const std::shared_ptr<FmEval_Chebyshev3>& instance(int m_max, double = 0.0) {

        // thread-safe per C++11 standard [6.7.4]
        static std::shared_ptr<FmEval_Chebyshev3> instance_ = 0;

        const bool need_new_instance = !instance_ || (instance_ && instance_->max_m() < m_max);
        if (need_new_instance) {
          auto new_instance = std::make_shared<FmEval_Chebyshev3>(m_max);
          instance_ = new_instance; // thread-safe
        }

        return instance_;
      }
#endif

      /// @return the maximum value of m for which the Boys function can be computed with this object
      int max_m() const { return mmax; }

      /// fills in Fm with computed Boys function values for m in [0,mmax]
      /// @param[out] Fm array to be filled in with the Boys function values, must be at least mmax+1 elements long
      /// @param[in] x the Boys function argument
      /// @param[in] mmax the maximum value of m for which Boys function will be computed; mmax must be <= the value returned by max_m
      inline void eval(Real* Fm, Real x, int m_max) const {

        // large T => use upward recursion
        // cost = 1 div + 1 sqrt + (1 + 2*(m-1)) muls
        if (x > T_crit) {
          const double one_over_x = 1.0/x;
          Fm[0] = 0.88622692545275801365 * sqrt(one_over_x); // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
          if (m_max == 0)
            return;
          // this upward recursion formula omits - e^(-x)/(2x), which for x>T_crit is <1e-15
          for (int i = 1; i <= m_max; i++)
            Fm[i] = Fm[i - 1] * numbers_.ihalf[i] * one_over_x; // see Eq. (9.8.13)
          return;
        }

        // ---------------------------------------------
        // small and intermediate arguments => interpolate Fm and (optional) downward recursion
        // ---------------------------------------------
        // about which point on the grid to interpolate?
        const double xd = x * one_over_delta;
        const int iv = int(xd); // the interval
        // INTERPOLATION_AND_RECURSION== true? evaluate by interpolation for LARGEST m only
        // INTERPOLATION_AND_RECURSION==false? evaluate by interpolation for ALL m
        const int m_min = INTERPOLATION_AND_RECURSION ? m_max : 0;

#if defined(__AVX__) || defined(__SSE2__)
        const auto x2 = xd*xd;
        const auto x3 = x2*xd;
#  if defined (__AVX__)
        libint2::simd::VectorAVXDouble xvec(1., xd, x2, x3);
#  else // defined(__SSE2__)
        libint2::simd::VectorSSEDouble x0vec(1., xd);
        libint2::simd::VectorSSEDouble x1vec(x2, x3);
#  endif
#endif // SSE2 || AVX

        const Real *d = c + INTERPOLATION_ORDER * (iv * (mmax+1) + m_min); // ptr to the interpolation data for m=mmin
        int m = m_min;
#if defined(__AVX__)
        if (m_max-m >=3) {
          const int unroll_size = 4;
          const int m_fence = (m_max + 2 - unroll_size);
          for(; m<m_fence; m+=unroll_size, d+=INTERPOLATION_ORDER*unroll_size) {
            libint2::simd::VectorAVXDouble d0v, d1v, d2v, d3v;
            d0v.load_aligned(d);
            d1v.load_aligned(d+INTERPOLATION_ORDER);
            d2v.load_aligned(d+2*INTERPOLATION_ORDER);
            d3v.load_aligned(d+3*INTERPOLATION_ORDER);
            libint2::simd::VectorAVXDouble fm0 = d0v * xvec;
            libint2::simd::VectorAVXDouble fm1 = d1v * xvec;
            libint2::simd::VectorAVXDouble fm2 = d2v * xvec;
            libint2::simd::VectorAVXDouble fm3 = d3v * xvec;
            libint2::simd::VectorAVXDouble sum0123 = horizontal_add(fm0, fm1, fm2, fm3);
            sum0123.convert(&Fm[m]);
          }
        } // unroll_size=4
        if (m_max-m >=1) {
          const int unroll_size = 2;
          const int m_fence = (m_max + 2 - unroll_size);
          for(; m<m_fence; m+=unroll_size, d+=INTERPOLATION_ORDER*unroll_size) {
            libint2::simd::VectorAVXDouble d0v, d1v;
            d0v.load_aligned(d);
            d1v.load_aligned(d+INTERPOLATION_ORDER);
            libint2::simd::VectorAVXDouble fm0 = d0v * xvec;
            libint2::simd::VectorAVXDouble fm1 = d1v * xvec;
            libint2::simd::VectorSSEDouble sum01 = horizontal_add(fm0, fm1);
            sum01.convert(&Fm[m]);
          }
        } // unroll_size=2
        { // no unrolling
          for(; m<=m_max; ++m, d+=INTERPOLATION_ORDER) {
            libint2::simd::VectorAVXDouble dvec;
            dvec.load_aligned(d);
            libint2::simd::VectorAVXDouble fm_prereduce = dvec * xvec;
            Fm[m] = horizontal_add(fm_prereduce);
          }
        }
#elif defined(__SSE2__)
        if (m_max-m >=1) {
          const int unroll_size = 2;
          const int m_fence = (m_max + 2 - unroll_size);
          for(; m<m_fence; m+=unroll_size, d+=INTERPOLATION_ORDER*unroll_size) {
            libint2::simd::VectorSSEDouble d00v, d01v, d10v, d11v;
            d00v.load_aligned(d);
            d01v.load_aligned(d+2);
            d10v.load_aligned(d+4); // d + INTERPOLATION_ORDER
            d11v.load_aligned(d+6);
            libint2::simd::VectorSSEDouble fm00 = d00v * x0vec;
            libint2::simd::VectorSSEDouble fm01 = d01v * x1vec;
            libint2::simd::VectorSSEDouble fm10 = d10v * x0vec;
            libint2::simd::VectorSSEDouble fm11 = d11v * x1vec;
            libint2::simd::VectorSSEDouble sum01 = horizontal_add(fm00, fm10) + horizontal_add(fm01, fm11);
            sum01.convert(&Fm[m]);
          }
        } // unroll_size=2
        { // no unrolling
          for(; m<=m_max; ++m, d+=INTERPOLATION_ORDER) {
            libint2::simd::VectorSSEDouble d0vec, d1vec;
            d0vec.load_aligned(d);
            d1vec.load_aligned(d+2);
            Fm[m] = horizontal_add(d0vec * x0vec + d1vec * x1vec);
          }
        }
#else // not SSE2 nor AVX available
        for(int m=m_min; m<=m_max; ++m, d+=INTERPOLATION_ORDER) {
          Fm[m] = d[0]
                + xd * (d[1] + xd * (d[2] + xd * d[3]));

          //        // check against the reference value
          //        if (false) {
          //          double refvalue = FmEval_Reference2<double>::eval(x, m, 1e-15); // compute F(T)
          //          if (abs(refvalue - Fm[m]) > 1e-10) {
          //            std::cout << "T = " << x << " m = " << m << " cheb = "
          //                      << Fm[m] << " ref = " << refvalue << std::endl;
          //          }
          //        }
        }
#endif


        // use downward recursion (Eq. (9.8.14) in HJO)
        // WARNING: do not turn this on ... on modern CPU exp is very slow
        if (INTERPOLATION_AND_RECURSION && m_max > 0) {
          const bool INTERPOLATION_AND_RECURSION_is_slow_on_modern_CPU = false;
          assert(INTERPOLATION_AND_RECURSION_is_slow_on_modern_CPU);
          const Real x2 = 2.0 * x;
          const Real exp_x = exp(-x);
          for (int m = m_max - 1; m >= 0; m--)
            Fm[m] = (Fm[m + 1] * x2 + exp_x) * numbers_.twoi1[m];
        }
      }

    private:

      /* ----------------------------------------------------------------------------
       This function here creates the expansion coefficients for a single interval

       ON INPUT  a,b  : the interval boundaries
       cc   : a pointer to the appropriate place in
       the coefficient table
       m    : the F[m] to generate
       ON OUTPUT cc   : cc[0]-cc[3] hold the coefficients
       ---------------------------------------------------------------------------- */
      void MakeCoeffs(double a, double b, Real *cc, int m) {
        int k, j;
        double f[128], ac[128], Fm[128];
        double sum;
        // characterize the interval
        double TwoDelta = b - a;
        double Delta = 0.5 * TwoDelta;
        double HalfDelta = 0.5 * Delta;
        double XXX = a + Delta;

        const double absolute_precision = 1e-100; // compute as precisely as possible
        FmEval_Reference2<double>::eval(Fm, XXX, m + INTERPOLATION_ORDER + 20,
                                        absolute_precision);

        for (k = 0; k <= INTERPOLATION_ORDER + 20; k++) {
          if ((k % 2) == 0)
            f[k] = Fm[k + m];
          else
            f[k] = -Fm[k + m];
        }
        // calculate the coefficients a
        double fac;
        for (j = 0; j < INTERPOLATION_ORDER; j++) {
          if (j == 0)
            fac = 1.0;
          else
            fac = 2.0 * pow(HalfDelta, (double) j);
          sum = 0.0;
          for (k = 0; k < 10; k++)
            sum += f[j + 2 * k] * pow(HalfDelta, (double) (2 * k)) / numbers_.fac[k]
                / numbers_.fac[k + j];
          ac[j] = fac * sum;
        }
        // calculate the coefficients c that are Gill's f's
        double arg = -XXX / Delta;
        double arg2 = arg * arg;
        double arg3 = arg2 * arg;
        auto cc0 = (ac[0] - ac[2]) + (ac[1] - 3.0 * ac[3]) * arg
            + 2.0 * ac[2] * arg2 + 4.0 * ac[3] * arg3;

        auto cc1 = (2.0 * ac[1] - 6.0 * ac[3]) + 8.0 * ac[2] * arg
            + 24.0 * ac[3] * arg2;
        auto cc2 = 8.0 * ac[2] + 48.0 * ac[3] * arg;
        auto cc3 = 32.0 * ac[3];
        cc[0] = cc0;
        cc[1] = cc1;
        cc[2] = cc2;
        cc[3] = cc3;
      }

      /* ----------------------------------------------------------------------------
       This function makes the expansion coefficients for all intervals


       ON INPUT  m    : the highest F[m] to generate

       ON OUTPUT  c   : the coefficients c[m][i] are generated
       ---------------------------------------------------------------------------- */
      void init() {
        int iv, im;

        // get memory
        void* result;
        posix_memalign(&result, 4*sizeof(Real), (mmax + 1) * NGRID * INTERPOLATION_ORDER * sizeof(Real));
        c = static_cast<Real*>(result);

        // make expansion coefficients for each grid value of T
        for (iv = 0; iv < NGRID; iv++) {
          const auto a = iv * delta;
          const auto b = a + delta;

          // loop over all m values and make the coefficients
          for (im = 0; im <= mmax; im++) {
            MakeCoeffs(a, b, c + (iv * (mmax+1) + im) * INTERPOLATION_ORDER, im);
          }
        }
      }

  };

#ifndef STATIC_OON
#define STATIC_OON
  namespace {
    const double oon[] = {0.0, 1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0, 1.0/5.0, 1.0/6.0, 1.0/7.0, 1.0/8.0, 1.0/9.0, 1.0/10.0, 1.0/11.0};
  }
#endif

  /** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \, {\rm d}u \$,
    * using Taylor interpolation of up to 8-th order.
    * @tparam Real the type to use for all floating-point computations. Must support std::exp, std::pow, std::fabs, std::max, and std::floor.
    * @tparam INTERPOLATION_ORDER the interpolation order. The higher the order the less memory this object will need, but the computational cost will increase (usually very slightly)
    */
  template<typename Real = double, int INTERPOLATION_ORDER = 7>
  class FmEval_Taylor {
    public:
      static const int max_interp_order = 8;
      static const bool INTERPOLATION_AND_RECURSION = false; // compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only
      const double soft_zero_;

      /// Constructs the object to be able to compute Boys funcion for m in [0,mmax], with relative \c precision
      FmEval_Taylor(unsigned int mmax, Real precision) :
          soft_zero_(1e-6), cutoff_(precision), numbers_(
              INTERPOLATION_ORDER + 1, 2 * (mmax + INTERPOLATION_ORDER - 1)) {

        assert(mmax <= 63);

        const double sqrt_pi = std::sqrt(M_PI);

        /*---------------------------------------
         We are doing Taylor interpolation with
         n=TAYLOR_ORDER terms here:
         error <= delT^n/(n+1)!
         ---------------------------------------*/
        delT_ = 2.0
            * std::pow(cutoff_ * numbers_.fac[INTERPOLATION_ORDER + 1],
                       1.0 / INTERPOLATION_ORDER);
        oodelT_ = 1.0 / delT_;
        max_m_ = mmax + INTERPOLATION_ORDER - 1;

        T_crit_ = new Real[max_m_ + 1]; /*--- m=0 is included! ---*/
        max_T_ = 0;
        /*--- Figure out T_crit for each m and put into the T_crit ---*/
        for (int m = max_m_; m >= 0; --m) {
          /*------------------------------------------
           Damped Newton-Raphson method to solve
           T^{m-0.5}*exp(-T) = epsilon*Gamma(m+0.5)
           The solution is the max T for which to do
           the interpolation
           ------------------------------------------*/
          double T = -log(cutoff_);
          const double egamma = cutoff_ * sqrt_pi * numbers_.df[2 * m]
              / std::pow(2.0, m);
          double T_new = T;
          double func;
          do {
            const double damping_factor = 0.2;
            T = T_new;
            /* f(T) = the difference between LHS and RHS of the equation above */
            func = std::pow(T, m - 0.5) * std::exp(-T) - egamma;
            const double dfuncdT = ((m - 0.5) * std::pow(T, m - 1.5)
                - std::pow(T, m - 0.5)) * std::exp(-T);
            /* f(T) has 2 roots and has a maximum in between. If f'(T) > 0 we are to the left of the hump. Make a big step to the right. */
            if (dfuncdT > 0.0) {
              T_new *= 2.0;
            } else {
              /* damp the step */
              double deltaT = -func / dfuncdT;
              const double sign_deltaT = (deltaT > 0.0) ? 1.0 : -1.0;
              const double max_deltaT = damping_factor * T;
              if (std::fabs(deltaT) > max_deltaT)
                deltaT = sign_deltaT * max_deltaT;
              T_new = T + deltaT;
            }
            if (T_new <= 0.0) {
              T_new = T / 2.0;
            }
          } while (std::fabs(func / egamma) >= soft_zero_);
          T_crit_[m] = T_new;
          const int T_idx = (int) std::floor(T_new / delT_);
          max_T_ = std::max(max_T_, T_idx);
        }

        // allocate the grid (see the comments below)
        {
          const int nrow = max_T_ + 1;
          const int ncol = max_m_ + 1;
          grid_ = new Real*[nrow];
          grid_[0] = new Real[nrow * ncol];
          //std::cout << "Allocated interpolation table of " << nrow * ncol << " reals" << std::endl;
          for (int r = 1; r < nrow; ++r)
            grid_[r] = grid_[r - 1] + ncol;
        }

        /*-------------------------------------------------------
         Tabulate the gamma function from t=delT to T_crit[m]:
         1) include T=0 though the table is empty for T=0 since
         Fm(0) is simple to compute
         -------------------------------------------------------*/
        /*--- do the mmax first ---*/
        for (int T_idx = max_T_; T_idx >= 0; --T_idx) {
          const double T = T_idx * delT_;
          libint2::FmEval_Reference2<double>::eval(grid_[T_idx], T, max_m_, 1e-100);
        }
      }

      ~FmEval_Taylor() {
        delete[] T_crit_;
        delete[] grid_[0];
        delete[] grid_;
      }

      // some features require at least C++11
#if __cplusplus > 199711L
      /// Singleton interface allows to manage the lone instance;
      /// adjusts max m and precision values as needed in thread-safe fashion
      static const std::shared_ptr<FmEval_Taylor>& instance(unsigned int mmax, Real precision) {

        // thread-safe per C++11 standard [6.7.4]
        static std::shared_ptr<FmEval_Taylor> instance_ = 0;

        const bool need_new_instance = !instance_ ||
                                       (instance_ && (instance_->max_m() < mmax ||
                                                      instance_->precision() > precision));
        if (need_new_instance) {
          auto new_instance = std::make_shared<FmEval_Taylor>(mmax, precision);
          instance_ = new_instance; // thread-safe
        }

        return instance_;
      }
#endif

      /// @return the maximum value of m for which this object can compute the Boys function
      int max_m() const { return max_m_ - INTERPOLATION_ORDER + 1; }
      /// @return the precision with which this object can compute the Boys function
      Real precision() const { return cutoff_; }

      /// computes Boys function values with m index in range [0,mmax]
      /// @param[out] Fm array to be filled in with the Boys function values, must be at least mmax+1 elements long
      /// @param[in] x the Boys function argument
      /// @param[in] mmax the maximum value of m for which Boys function will be computed;
      ///                  it must be <= the value returned by max_m() (this is not checked)
      void eval(Real* Fm, Real T, int mmax) const {
        const double sqrt_pio2 = 1.2533141373155002512;
        const double two_T = 2.0 * T;

        // stop recursion at mmin
        const int mmin = INTERPOLATION_AND_RECURSION ? mmax : 0;
        /*-------------------------------------
         Compute Fm(T) from mmax down to mmin
         -------------------------------------*/
        const bool use_upward_recursion = true;
        if (use_upward_recursion) {
//          if (T > 30.0) {
          if (T > T_crit_[0]) {
            const double one_over_x = 1.0/T;
            Fm[0] = 0.88622692545275801365 * sqrt(one_over_x); // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
            if (mmax == 0)
              return;
            // this upward recursion formula omits - e^(-x)/(2x), which for x>T_crit is <1e-15
            for (int i = 1; i <= mmax; i++)
              Fm[i] = Fm[i - 1] * numbers_.ihalf[i] * one_over_x; // see Eq. (9.8.13)
            return;
          }
        }

        // since Tcrit grows with mmax, this condition only needs to be determined once
        if (T > T_crit_[mmax]) {
          double pow_two_T_to_minusjp05 = std::pow(two_T, -mmax - 0.5);
          for (int m = mmax; m >= mmin; --m) {
            /*--- Asymptotic formula ---*/
            Fm[m] = numbers_.df[2 * m] * sqrt_pio2 * pow_two_T_to_minusjp05;
            pow_two_T_to_minusjp05 *= two_T;
          }
        }
        else
        {
          const int T_ind = (int) (0.5 + T * oodelT_);
          const Real h = T_ind * delT_ - T;
          const Real* F_row = grid_[T_ind] + mmin;

#if defined (__AVX__)
          libint2::simd::VectorAVXDouble h0123, h4567;
          if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
            const double h2 = h*h*oon[2];
            const double h3 = h2*h*oon[3];
            h0123 = libint2::simd::VectorAVXDouble (1.0, h, h2, h3);
            if (INTERPOLATION_ORDER == 7) {
              const double h4 = h3*h*oon[4];
              const double h5 = h4*h*oon[5];
              const double h6 = h5*h*oon[6];
              const double h7 = h6*h*oon[7];
              h4567 = libint2::simd::VectorAVXDouble (h4, h5, h6, h7);
            }
          }
          //          libint2::simd::VectorAVXDouble h0123(1.0);
          //          libint2::simd::VectorAVXDouble h4567(1.0);
#endif

          int m = mmin;
          if (mmax-m >=1) {
            const int unroll_size = 2;
            const int m_fence = (mmax + 2 - unroll_size);
            for(; m<m_fence; m+=unroll_size, F_row+=unroll_size) {

#if defined(__AVX__)
              if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
                 libint2::simd::VectorAVXDouble fr0_0123; fr0_0123.load(F_row);
                 libint2::simd::VectorAVXDouble fr1_0123; fr1_0123.load(F_row+1);
                 libint2::simd::VectorSSEDouble fm01 = horizontal_add(fr0_0123*h0123, fr1_0123*h0123);
                 if (INTERPOLATION_ORDER == 7) {
                   libint2::simd::VectorAVXDouble fr0_4567; fr0_4567.load(F_row+4);
                   libint2::simd::VectorAVXDouble fr1_4567; fr1_4567.load(F_row+5);
                   fm01 += horizontal_add(fr0_4567*h4567, fr1_4567*h4567);
                 }
                 fm01.convert(&Fm[m]);
              }
              else {
#endif
              Real total0 = 0.0, total1 = 0.0;
              for(int i=INTERPOLATION_ORDER; i>=1; --i) {
                total0 = oon[i]*h*(F_row[i] + total0);
                total1 = oon[i]*h*(F_row[i+1] + total1);
              }
              Fm[m] = F_row[0] + total0;
              Fm[m+1] = F_row[1] + total1;
#if defined(__AVX__)
              }
#endif
            }
          } // unroll_size = 2
          if (m<=mmax) { // unroll_size = 1
#if defined(__AVX__)
            if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
              libint2::simd::VectorAVXDouble fr0123; fr0123.load(F_row);
              if (INTERPOLATION_ORDER == 7) {
                libint2::simd::VectorAVXDouble fr4567; fr4567.load(F_row+4);
                libint2::simd::VectorSSEDouble fm = horizontal_add(fr0123*h0123, fr4567*h4567);
                Fm[m] = horizontal_add(fm);
              }
              else { // INTERPOLATION_ORDER == 3
                Fm[m] = horizontal_add(fr0123*h0123);
              }
            }
            else {
#endif
            Real total = 0.0;
            for(int i=INTERPOLATION_ORDER; i>=1; --i) {
              total = oon[i]*h*(F_row[i] + total);
            }
            Fm[m] = F_row[0] + total;
#if defined(__AVX__)
            }
#endif
          } // unroll_size = 1

          // check against the reference value
//          if (false) {
//            double refvalue = FmEval_Reference2<double>::eval(T, mmax, 1e-15); // compute F(T) with m=mmax
//            if (abs(refvalue - Fm[mmax]) > 1e-14) {
//              std::cout << "T = " << T << " m = " << mmax << " cheb = "
//                  << Fm[mmax] << " ref = " << refvalue << std::endl;
//            }
//          }

        } // if T < T_crit

        /*------------------------------------
         And then do downward recursion in j
         ------------------------------------*/
        if (INTERPOLATION_AND_RECURSION && mmin > 0) {
          const Real exp_mT = std::exp(-T);
          for (int m = mmin - 1; m >= 0; --m) {
            Fm[m] = (exp_mT + two_T * Fm[m+1]) * numbers_.twoi1[m];
          }
        }
      }

    private:
      Real **grid_; /* Table of "exact" Fm(T) values. Row index corresponds to
       values of T (max_T+1 rows), column index to values
       of m (max_m+1 columns) */
      Real delT_; /* The step size for T, depends on cutoff */
      Real oodelT_; /* 1.0 / delT_, see above */
      Real cutoff_; /* Tolerance cutoff used in all computations of Fm(T) */
      int max_m_; /* Maximum value of m in the table, depends on cutoff
       and the number of terms in Taylor interpolation */
      int max_T_; /* Maximum index of T in the table, depends on cutoff
       and m */
      Real *T_crit_; /* Maximum T for each row, depends on cutoff;
       for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
       for a given m and T_idx > max_T_idx[m] use the asymptotic formula */

      ExpensiveNumbers<double> numbers_;

      /**
       * Power series estimate of the error introduced by replacing
       * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with analytically
       * integrable \f$ \int_0^\infty \exp(-T t^2) t^{2 m} \, \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
       * @param m
       * @param T
       * @return the error estimate
       */
      static double truncation_error(unsigned int m, double T) {
        const double m2= m * m;
        const double m3= m2 * m;
        const double m4= m2 * m2;
        const double T2= T * T;
        const double T3= T2 * T;
        const double T4= T2 * T2;
        const double T5= T2 * T3;

        const double result = exp(-T) * (105 + 16*m4 + 16*m3*(T - 8) - 30*T + 12*T2
            - 8*T3 + 16*T4 + 8*m2*(43 - 9*T + 2*T2) +
            4*m*(-88 + 23*T - 8*T2 + 4*T3))/(32*T5);
        return result;
      }
      /**
       * Leading-order estimate of the error introduced by replacing
       * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with analytically
       * integrable \f$ \int_0^\infty \exp(-T t^2) t^{2 m} \, \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
       * @param m
       * @param T
       * @return the error estimate
       */
      static double truncation_error(double T) {
        const double result = exp(-T) /(2*T);
        return result;
      }
  };


  //////////////////////////////////////////////////////////
  /// core integral for Yukawa and exponential interactions
  //////////////////////////////////////////////////////////

#if 0
  /**
   * Evaluates core integral for the Yukawa potential \f$ \exp(- \zeta r) / r \f$
   * @tparam Real real type
   */
  template<typename Real>
  struct YukawaGmEval {

      static const int mmin = -1;

      ///
      YukawaGmEval(unsigned int mmax, Real precision) :
        mmax_(mmax), precision_(precision),
        numbers_(),
        Gm_0_U_(256) // should be enough to hold up to G_{255}(0,U)
      { }

      unsigned int max_m() const { return mmax; }
      /// @return the precision with which this object can compute the result
      Real precision() const { return precision_; }

      ///
      void eval_yukawa(Real* Gm, Real T, Real U, size_t mmax, Real absolute_precision) {
        assert(false); // not yet implemented
      }
      ///
      void eval_slater(Real* Gm, Real T, Real U, size_t mmax, Real absolute_precision) {
        assert(false); // not yet implemented
      }

      /// Scheme 1 of Ten-no: upward recursion from \f$ G_{-1} (T,U) \f$ and \f$ G_0 (T,U) \f$
      /// T must be non-zero!
      /// @param[out] Gm \f$ G_m(T,U), m=-1..mmax \f$
      static void eval_yukawa_s1(Real* Gm, Real T, Real U, size_t mmax) {
        Real G_m1;

        const Real sqrt_U = sqrt(U);
        const Real sqrt_T = sqrt(T);
        const Real oo_sqrt_T = 1 / sqrt_T;
        const Real oo_sqrt_U = 1 / sqrt_U;
        const Real exp_mT = exp(-T);
        const Real kappa = sqrt_U - sqrt_T;
        const Real lambda = sqrt_U + sqrt_T;
        const Real sqrtPi_over_4(0.44311346272637900682454187083528629569938736403060);
        const Real pfac = sqrtPi_over_4 * exp_mT;
        const Real erfc_k = exp(kappa*kappa) * (1 - erf(kappa));
        const Real erfc_l = exp(lambda*lambda) * (1 - erf(lambda));

        Gm[0] = pfac * (erfc_k + erfc_l) * oo_sqrt_U;
        Gm[1] = pfac * (erfc_k - erfc_l) * oo_sqrt_T;
        if (mmax > 0) {

          // first application of URR
          const Real oo_two_T = 0.5 / T;
          const Real two_U = 2.0 * U;

          for(unsigned int m=1, two_m_minus_1=1; m<=mmax; ++m, two_m_minus_1+=2) {
            Gm[m+1] = oo_two_T * ( two_m_minus_1 * Gm[m] + two_U * Gm[m-1] - exp_mT);
          }
        }

        return;
      }

      /// Scheme 2 of Ten-no:
      /// - evaluate G_m(0,U) for m = mmax ... mmax+n, where n is the number of terms in Maclaurin expansion
      ///   how? see eval_yukawa_Gm0U
      /// - then MacLaurin expansion for \f$ G_{m_{\rm max}}(T,U) \f$ and \f$ G_{m_{\rm max}-1}(T,U) \f$
      /// - then downward recursion
      /// @param[out] Gm \f$ G_m(T,U), m=-1..mmax \f$
      void eval_yukawa_s2(Real* Gm, Real T, Real U, size_t mmax) {

        // TODO estimate the number of expansion terms for the given precision
        const int expansion_order = 60;
        eval_yukawa_Gm0U(Gm_0_U_, U, mmax - 1 + expansion_order);

        // Maclaurin


        // downward recursion
        //Gm[m + 1] = 1/(2 U) (E^-T - (2 m + 3) Gm[[m + 2]] + 2 T Gm[[m + 3]])
        const Real one_over_twoU = 0.5 / U;
        const Real one_over_twoU = 2.0 * T;
        const Real exp_mT = exp(-T);
        for(int m=mmax-2; m>=-1; --m)
          Gm[m] = one_over_twoU (exp_mT - numbers_.twoi1[m+1] * Gm[m+1] + twoT Gm[m+2])

        // testing ...
        std::copy(Gm_0_U_.begin()+1, Gm_0_U_.begin()+mmax+2, Gm);

        return;
      }

      /// Scheme 3 of Ten-no:
      /// - evaluate G_m(0,U) for m = 0 ... mmax+n, where n is the max order of terms in Maclaurin expansion
      ///   how? see eval_yukawa_Gm0U
      /// - then MacLaurin expansion for \f$ G_{m}(T,U) \f$ for m = 0 ... mmax
      /// @param[out] Gm \f$ G_m(T,U), m=-1..mmax \f$
      void eval_yukawa_s3(Real* Gm, Real T, Real U, size_t mmax) {

        // Ten-no's prescription:
        //

        assert(false);

        // testing ...
        std::copy(Gm_0_U_.begin()+1, Gm_0_U_.begin()+mmax+2, Gm);

        return;
      }


      /**
       * computes prerequisites for MacLaurin expansion of Gm(T,U)
       * for m in [-1,mmax); uses Ten-no's prescription, i.e.
       *
       *
       * @param[out] Gm0U
       * @param[in] U
       * @param[in] mmax
       */
      void eval_yukawa_Gm0U(Real* Gm0U, Real U, int mmax, int mmin = -1) {

        // Ten-no's prescription:
        // start with Gm*(0,T)
        // 1) for U < 5, m* = -1
        // 2) for U > 5, m* = min(U,mmax)
        int mstar;

        // G_{-1} (0,U) is easy
        if (U < 5.0) {
          mstar = -1;

          const Real sqrt_U = sqrt(U);
          const Real exp_U = exp(U);
          const Real oo_sqrt_U = 1 / sqrt_U;
          const Real sqrtPi_over_2(
              0.88622692545275801364908374167057259139877472806119);
          const Real pfac = sqrtPi_over_2 * exp_U;
          const Real erfc_sqrt_U = 1.0 - erf(sqrt_U);
          Gm_0_U_[0] = pfac * exp_U * oo_sqrt_U * erfc_sqrt_U;
          // can get G0 for "free"
          // this is the l'Hopital-transformed expression for G_0 (0,T)
//          const Real sqrtPi(
//              1.7724538509055160272981674833411451827975494561224);
//          Gm_0_U_[1] = 1.0 - exp_U * sqrtPi * sqrt_U * erfc_sqrt_U;
        }
        else { // use continued fraction for m*
          mstar = std::min((size_t)U,(size_t)mmax);
          const bool implemented = false;
          assert(implemented == true);
        }

        { // use recursion if needed
          const Real two_U = 2.0 * U;
          // simplified URR
          if (mmax > mstar) {
            for(int m=mstar+1; m<=mmax; ++m) {
              Gm_0_U_[m+1] = numbers_.twoi1[m] * (1.0 - two_U * Gm_0_U_[m]);
            }
          }

          // simplified DRR
          if (mstar > mmin) { // instead of -1 because we trigger this only for U > 5
            const Real one_over_U = 2.0 / two_U;
            for(int m=mstar-1; m>=mmin; --m) {
              Gm_0_U_[m+1] = one_over_U * ( 0.5 - numbers_.ihalf[m+2] * Gm_0_U_[m+2]);
            }
          }
        }

        // testing ...
        std::copy(Gm_0_U_.begin()+1, Gm_0_U_.begin()+mmax+2, Gm0U);

        return;
      }


      /// computes a single value of G_{-1}(T,U)
      static Real eval_Gm1(Real T, Real U) {
        const Real sqrt_U = sqrt(U);
        const Real sqrt_T = sqrt(T);
        const Real exp_mT = exp(-T);
        const Real kappa = sqrt_U - sqrt_T;
        const Real lambda = sqrt_U + sqrt_T;
        const Real sqrtPi_over_4(0.44311346272637900682454187083528629569938736403060);
        const Real result = sqrtPi_over_4 * exp_mT *
            (exp(kappa*kappa) * (1 - erf(kappa)) + exp(lambda*lambda) * (1 - erf(lambda))) / sqrt_U;
        return result;
      }
      /// computes a single value of G_0(T,U)
      static Real eval_G0(Real T, Real U) {
        const Real sqrt_U = sqrt(U);
        const Real sqrt_T = sqrt(T);
        const Real exp_mT = exp(-T);
        const Real kappa = sqrt_U - sqrt_T;
        const Real lambda = sqrt_U + sqrt_T;
        const Real sqrtPi_over_4(0.44311346272637900682454187083528629569938736403060);
        const Real result = sqrtPi_over_4 * exp_mT *
            (exp(kappa*kappa) * (1 - erf(kappa)) - exp(lambda*lambda) * (1 - erf(lambda))) / sqrt_T;
        return result;
      }
      /// computes \f$ G_{-1}(T,U) \f$ and \f$ G_{0}(T,U) \f$ , both are needed for Yukawa and Slater integrals
      /// @param[out] result result[0] contains \f$ G_{-1}(T,U) \f$, result[1] contains \f$ G_{0}(T,U) \f$
      static void eval_G_m1_0(Real* result, Real T, Real U) {
        const Real sqrt_U = sqrt(U);
        const Real sqrt_T = sqrt(T);
        const Real oo_sqrt_U = 1 / sqrt_U;
        const Real oo_sqrt_T = 1 / sqrt_T;
        const Real exp_mT = exp(-T);
        const Real kappa = sqrt_U - sqrt_T;
        const Real lambda = sqrt_U + sqrt_T;
        const Real sqrtPi_over_4(0.44311346272637900682454187083528629569938736403060);
        const Real pfac = sqrtPi_over_4 * exp_mT;
        const Real erfc_k = exp(kappa*kappa) * (1 - erf(kappa));
        const Real erfc_l = exp(lambda*lambda) * (1 - erf(lambda));
        result[0] = pfac * (erfc_k + erfc_l) * oo_sqrt_U;
        result[1] = pfac * (erfc_k - erfc_l) * oo_sqrt_T;
      }

      /// computes a single value of G(T,U) using MacLaurin series.
      static Real eval_MacLaurinT(Real T, Real U, size_t m, Real absolute_precision) {
        assert(false); // not yet implemented
        return 0.0;
      }

    private:
      std::vector<Real> Gm_0_U_; // used for MacLaurin expansion
      unsigned int mmax_;
      Real precision_;
      ExpensiveNumbers<Real> numbers_;

      // since evaluation may involve several functions, will store some intermediate constants here
      // to avoid the cost of extra parameters
      //Real exp_U_;
      //Real exp_mT_;

      size_t count_tenno_algorithm_branches[3]; // counts the number of times each branch Ten-no algorithm
                                                // was picked

  };
#endif

  //////////////////////////////////////////////////////////
  /// core integrals r12^k \sum_i \exp(- a_i r_12^2)
  //////////////////////////////////////////////////////////

  /// each thread needs its own scratch if k==-1
  template <typename Real, int k> struct GaussianGmEvalScratch;
  /// create this object for each thread that needs to use GaussianGmEval<Real,-1>
  template <typename Real>
  struct GaussianGmEvalScratch<Real, -1> {
      std::vector<Real> Fm_;
      std::vector<Real> g_i;
      std::vector<Real> r_i;
      std::vector<Real> oorhog_i;
      GaussianGmEvalScratch() {}
      GaussianGmEvalScratch(int mmax) {
        init(mmax);
      }
      void init(int mmax) {
        Fm_.resize(mmax+1);
        g_i.resize(mmax+1);
        r_i.resize(mmax+1);
        oorhog_i.resize(mmax+1);
        g_i[0] = 1.0;
        r_i[0] = 1.0;
      }
  };


  /**
   * Evaluates core integral \$ G_m(\rho, T) = \left( - \frac{\partial}{\partial T} \right)^n G_0(\rho,T) \f$,
   * \f$ G_0(\rho,T) = \int \exp(-\rho |\vec{r} - \vec{P} + \vec{Q}|^2) g(r) \, {\rm d}\vec{r} \f$
   * over a general contracted
   * Gaussian geminal \f$ g(r_{12}) = r_{12}^k \sum_i c_i \exp(- a_i r_{12}^2), \quad k = -1, 0, 2 \f$ .
   * The integrals are needed in R12/F12 methods with STG-nG correlation factors.
   * Specifically, for a correlation factor \f$ f(r_{12}) = \sum_i c_i \exp(- a_i r_{12}^2) \f$
   * integrals with the following kernels are needed:
   * <ul>
   *   <li> \f$ f(r_{12}) \f$  (k=0) </li>
   *   <li> \f$ f(r_{12}) / r_{12} \f$  (k=-1) </li>
   *   <li> \f$ f(r_{12})^2 \f$ (k=0, @sa GaussianGmEval::eval ) </li>
   *   <li> \f$ [f(r_{12}), [\hat{T}_1, f(r_{12})]] \f$ (k=2, @sa GaussianGmEval::eval ) </li>
   * </ul>
   *
   * N.B. ``Asymmetric'' kernels, \f$ f(r_{12}) g(r_{12}) \f$ and
   *   \f$ [f(r_{12}), [\hat{T}_1, g(r_{12})]] \f$, where f and g are two different geminals,
   *   can also be handled straightforwardly.
   *
   * \note for more details see DOI: 10.1039/b605188j
   */
  template<typename Real, int k>
  struct GaussianGmEval {

      /**
       * @param[in] mmax the evaluator will be used to compute Gm(T) for 0 <= m <= mmax
       */
      GaussianGmEval(int mmax, Real precision) : mmax_(mmax),
          precision_(precision), fm_eval_(0),
          numbers_(-1,-1,mmax) {
        assert(k == -1 || k == 0 || k == 2);
        // for k=-1 need to evaluate the Boys function
        if (k == -1) {
          fm_eval_ = new FmEval_Taylor<Real>(mmax_, precision_);
//          fm_eval_ = new FmEval_Chebyshev3(mmax_);
          scratch_.init(mmax_);
        }
      }

      ~GaussianGmEval() {
        delete fm_eval_;
        fm_eval_ = 0;
      }

      // some features require at least C++11
#if __cplusplus > 199711L
      /// Singleton interface allows to manage the lone instance;
      /// adjusts max m and precision values as needed in thread-safe fashion
      static const std::shared_ptr<GaussianGmEval>& instance(unsigned int mmax, Real precision) {

        // thread-safe per C++11 standard [6.7.4]
        static std::shared_ptr<GaussianGmEval> instance_ = 0;

        const bool need_new_instance = !instance_ ||
                                       (instance_ && (instance_->max_m() < mmax ||
                                                      instance_->precision() > precision));
        if (need_new_instance) {
          auto new_instance = std::make_shared<GaussianGmEval>(mmax, precision);
          instance_ = new_instance; // thread-safe
        }

        return instance_;
      }
#endif

      /// @return the maximum value of m for which the \f$ G_m(\rho, T) \f$ can be computed with this object
      int max_m() const { return mmax_; }
      /// @return the precision with which this object can compute the Boys function
      Real precision() const { return precision_; }

      /** computes \f$ G_m(\rho, T) \f$ using downward recursion.
       *
       * @warning NOT reentrant if \c k == -1 and C++11 is not available
       *
       * @param[out] Gm array to be filled in with the \f$ Gm(\rho, T) \f$ values, must be at least mmax+1 elements long
       * @param[in] rho
       * @param[in] T
       * @param[in] mmax mmax the maximum value of m for which Boys function will be computed;
       *                 it must be <= the value returned by max_m() (this is not checked)
       * @param[in] geminal the Gaussian geminal for which the core integral \f$ Gm(\rho, T) \f$ is computed
       * @param[in] scr if \c k ==-1 and need this to be reentrant, must provide ptr to the per-thread \c GaussianGmEvalScratch<Real,-1> object;
       *                no need to specify \c scr otherwise
       */
      template <typename AnyReal>
      void eval(Real* Gm, Real rho, Real T, size_t mmax,
                const std::vector<std::pair<AnyReal, AnyReal> >& geminal,
                void* scr = 0) {

        std::fill(Gm, Gm+mmax+1, Real(0));

        const double sqrt_rho = sqrt(rho);
        const double oo_sqrt_rho = 1.0/sqrt_rho;
        if (k == -1) {
          GaussianGmEvalScratch<Real, -1>& scratch = (scr == 0) ? scratch_ : *(reinterpret_cast<GaussianGmEvalScratch<Real, -1>*>(scr));
          for(int i=1; i<=mmax; i++) {
            scratch.r_i[i] = scratch.r_i[i-1] * rho;
          }
        }

        typedef typename std::vector<std::pair<AnyReal, AnyReal> >::const_iterator citer;
        const citer gend = geminal.end();
        for(citer i=geminal.begin(); i!= gend; ++i) {

          const double gamma = i->first;
          const double gcoef = i->second;
          const double rhog = rho + gamma;
          const double oorhog = 1.0/rhog;

          const double gorg = gamma * oorhog;
          const double rorg = rho * oorhog;
          const double sqrt_rho_org = sqrt_rho * oorhog;
          const double sqrt_rhog = sqrt(rhog);
          const double sqrt_rorg = sqrt_rho_org * sqrt_rhog;

          /// (ss|g12|ss)
          const Real const_SQRTPI_2(0.88622692545275801364908374167057259139877472806119); /* sqrt(pi)/2 */
          const double SS_K0G12_SS = gcoef * oo_sqrt_rho * const_SQRTPI_2 * rorg * sqrt_rorg * exp(-gorg*T);

          if (k == -1) {
            GaussianGmEvalScratch<Real, -1>& scratch = (scr == 0) ? scratch_ : *(reinterpret_cast<GaussianGmEvalScratch<Real, -1>*>(scr));

            const double rorgT = rorg * T;
            fm_eval_->eval(&scratch.Fm_[0], rorgT, mmax);

#if 1
            const Real const_2_SQRTPI(1.12837916709551257389615890312154517);   /* 2/sqrt(pi)     */
            const Real pfac = const_2_SQRTPI * sqrt_rhog * SS_K0G12_SS;
            scratch.oorhog_i[0] = pfac;
            for(int i=1; i<=mmax; i++) {
              scratch.g_i[i] = scratch.g_i[i-1] * gamma;
              scratch.oorhog_i[i] = scratch.oorhog_i[i-1] * oorhog;
            }
            for(int m=0; m<=mmax; m++) {
              Real ssss = 0.0;
              Real* bcm = numbers_.bc[m];
              for(int n=0; n<=m; n++) {
                ssss += bcm[n] * scratch.r_i[n] * scratch.g_i[m-n] * scratch.Fm_[n];
              }
              Gm[m] += ssss * scratch.oorhog_i[m];
            }
#endif
          }

          if (k == 0) {
            double ss_oper_ss_m = SS_K0G12_SS;
            Gm[0] += ss_oper_ss_m;
            for(int m=1; m<=mmax; ++m) {
              ss_oper_ss_m *= gorg;
              Gm[m] += ss_oper_ss_m;
            }
          }

          if (k == 2) {

            /// (ss|g12*r12^2|ss)
            const double rorgT = rorg * T;
            const double SS_K2G12_SS_0 = (1.5 + rorgT) * (SS_K0G12_SS * oorhog);
            const double SS_K2G12_SS_m1 = rorg * (SS_K0G12_SS * oorhog);

            double SS_K2G12_SS_gorg_m = SS_K2G12_SS_0 ;
            double SS_K2G12_SS_gorg_m1 = SS_K2G12_SS_m1;
            Gm[0] += SS_K2G12_SS_gorg_m;
            for(int m=1; m<=mmax; ++m) {
              SS_K2G12_SS_gorg_m *= gorg;
              Gm[m] += SS_K2G12_SS_gorg_m - m * SS_K2G12_SS_gorg_m1;
              SS_K2G12_SS_gorg_m1 *= gorg;
            }
          }

        }

      }

    private:
      int mmax_;
      Real precision_; //< absolute precision
      FmEval_Taylor<Real>* fm_eval_;
//      FmEval_Chebyshev3* fm_eval_;

      GaussianGmEvalScratch <Real, -1> scratch_; // only used in serial if k==-1

      ExpensiveNumbers<Real> numbers_;
  };

  /*
   *  Slater geminal fitting is available only if have LAPACK
   */
#if HAVE_LAPACK
  /*
  f[x_] := - Exp[-\[Zeta] x] / \[Zeta];

  ff[cc_, aa_, x_] := Sum[cc[[i]]*Exp[-aa[[i]] x^2], {i, 1, n}];
  */
  template <typename Real>
  Real
  fstg(Real zeta,
       Real x) {
    return -std::exp(-zeta*x)/zeta;
  }

  template <typename Real>
  Real
  fngtg(const std::vector<Real>& cc,
        const std::vector<Real>& aa,
        Real x) {
    Real value = 0.0;
    const Real x2 = x * x;
    const unsigned int n = cc.size();
    for(unsigned int i=0; i<n; ++i)
      value += cc[i] * std::exp(- aa[i] * x2);
    return value;
  }

  // --- weighting functions ---
  // L2 error is weighted by ww(x)
  // hence error is weighted by sqrt(ww(x))
  template <typename Real>
  Real
  wwtewklopper(Real x) {
    const Real x2 = x * x;
    return x2 * std::exp(-2 * x2);
  }
  template <typename Real>
  Real
  wwcusp(Real x) {
    const Real x2 = x * x;
    const Real x6 = x2 * x2 * x2;
    return std::exp(-0.005 * x6);
  }
  // default is Tew-Klopper
  template <typename Real>
  Real
  ww(Real x) {
    //return wwtewklopper(x);
    return wwcusp(x);
  }

  template <typename Real>
  Real
  norm(const std::vector<Real>& vec) {
    Real value = 0.0;
    const unsigned int n = vec.size();
    for(unsigned int i=0; i<n; ++i)
      value += vec[i] * vec[i];
    return value;
  }

  template <typename Real>
  void LinearSolveDamped(const std::vector<Real>& A,
                         const std::vector<Real>& b,
                         Real lambda,
                         std::vector<Real>& x) {
    const size_t n = b.size();
    std::vector<Real> Acopy(A);
    for(size_t m=0; m<n; ++m) Acopy[m*n + m]  *= (1 + lambda);
    std::vector<Real> e(b);

    //int info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, 1, &Acopy[0], n, &ipiv[0], &e[0], n );
    {
      std::vector<int> ipiv(n);
      int n = b.size();
      int one = 1;
      int info;
      dgesv_(&n, &one, &Acopy[0], &n, &ipiv[0], &e[0], &n, &info);
      assert (info == 0);
    }

    x = e;
  }

  /**
   * computes a least-squares fit of \f$ -exp(-\zeta r_{12})/\zeta = \sum_{i=1}^n c_i exp(-a_i r_{12}^2) \f$
   * on \f$ r_{12} \in [0, x_{\rm max}] \f$ discretized to npts.
   * @param[in] n
   * @param[in] zeta
   * @param[out] geminal
   * @param[in] xmin
   * @param[in] xmax
   * @param[in] npts
   */
  template <typename Real>
  void stg_ng_fit(unsigned int n,
                 Real zeta,
                 std::vector< std::pair<Real, Real> >& geminal,
                 Real xmin = 0.0,
                 Real xmax = 10.0,
                 unsigned int npts = 1001) {

    // initial guess
    std::vector<Real> cc(n, 1.0); // coefficients
    std::vector<Real> aa(n); // exponents
    for(unsigned int i=0; i<n; ++i)
      aa[i] = std::pow(3.0, (i + 2 - (n + 1)/2.0));

    // first rescale cc for ff[x] to match the norm of f[x]
    Real ffnormfac = 0.0;
    for(unsigned int i=0; i<n; ++i)
      for(unsigned int j=0; j<n; ++j)
        ffnormfac += cc[i] * cc[j]/std::sqrt(aa[i] + aa[j]);
    const Real Nf = std::sqrt(2.0 * zeta) * zeta;
    const Real Nff = std::sqrt(2.0) / (std::sqrt(ffnormfac) *
        std::sqrt(std::sqrt(M_PI)));
    for(unsigned int i=0; i<n; ++i) cc[i] *= -Nff/Nf;

    Real lambda0 = 1000; // damping factor is initially set to 1000, eventually should end up at 0
    const Real nu = 3.0; // increase/decrease the damping factor scale it by this
    const Real epsilon = 1e-15; // convergence
    const unsigned int maxniter = 200;

    // grid points on which we will fit
    std::vector<Real> xi(npts);
    for(unsigned int i=0; i<npts; ++i) xi[i] = xmin + (xmax - xmin)*i/(npts - 1);

    std::vector<Real> err(npts);

    const size_t nparams = 2*n; // params = expansion coefficients + gaussian exponents
    std::vector<Real> J( npts * nparams );
    std::vector<Real> delta(nparams);

//    std::cout << "iteration 0" << std::endl;
//    for(unsigned int i=0; i<n; ++i)
//      std::cout << cc[i] << " " << aa[i] << std::endl;

    Real errnormI;
    Real errnormIm1 = 1e3;
    bool converged = false;
    unsigned int iter = 0;
    while (!converged && iter < maxniter) {
//      std::cout << "Iteration " << ++iter << ": lambda = " << lambda0/nu << std::endl;

        for(unsigned int i=0; i<npts; ++i) {
          const Real x = xi[i];
          err[i] = (fstg(zeta, x) - fngtg(cc, aa, x)) * std::sqrt(ww(x));
        }
        errnormI = norm(err)/std::sqrt((Real)npts);

//        std::cout << "|err|=" << errnormI << std::endl;
        converged = std::abs((errnormI - errnormIm1)/errnormIm1) <= epsilon;
        if (converged) break;
        errnormIm1 = errnormI;

        for(unsigned int i=0; i<npts; ++i) {
          const Real x2 = xi[i] * xi[i];
          const Real sqrt_ww_x = std::sqrt(ww(xi[i]));
          const unsigned int ioffset = i * nparams;
          for(unsigned int j=0; j<n; ++j)
            J[ioffset+j] = (std::exp(-aa[j] * x2)) * sqrt_ww_x;
          const unsigned int ioffsetn = ioffset+n;
          for(unsigned int j=0; j<n; ++j)
            J[ioffsetn+j] = -  sqrt_ww_x * x2 * cc[j] * std::exp(-aa[j] * x2);
        }

        std::vector<Real> A( nparams * nparams);
        for(size_t r=0, rc=0; r<nparams; ++r) {
          for(size_t c=0; c<nparams; ++c, ++rc) {
            double Arc = 0.0;
            for(size_t i=0, ir=r, ic=c; i<npts; ++i, ir+=nparams, ic+=nparams)
              Arc += J[ir] * J[ic];
            A[rc] = Arc;
          }
        }

        std::vector<Real> b( nparams );
        for(size_t r=0; r<nparams; ++r) {
          Real br = 0.0;
          for(size_t i=0, ir=r; i<npts; ++i, ir+=nparams)
            br += J[ir] * err[i];
          b[r] = br;
        }

        // try decreasing damping first
        // if not successful try increasing damping until it results in a decrease in the error
        lambda0 /= nu;
        for(int l=-1; l<1000; ++l) {

          LinearSolveDamped(A, b, lambda0, delta );

          std::vector<double> cc_0(cc); for(unsigned int i=0; i<n; ++i) cc_0[i] += delta[i];
          std::vector<double> aa_0(aa); for(unsigned int i=0; i<n; ++i) aa_0[i] += delta[i+n];

          // if any of the exponents are negative the step is too large and need to increase damping
          bool step_too_large = false;
          for(unsigned int i=0; i<n; ++i)
            if (aa_0[i] < 0.0) {
              step_too_large = true;
              break;
            }
          if (!step_too_large) {
            std::vector<double> err_0(npts);
            for(unsigned int i=0; i<npts; ++i) {
              const double x = xi[i];
              err_0[i] = (fstg(zeta, x) - fngtg(cc_0, aa_0, x)) * std::sqrt(ww(x));
            }
            const double errnorm_0 = norm(err_0)/std::sqrt((double)npts);
            if (errnorm_0 < errnormI) {
              cc = cc_0;
              aa = aa_0;
              break;
            }
            else // step lead to increase of the error -- try dampening a bit more
              lambda0 *= nu;
          }
          else // too large of a step
            lambda0 *= nu;
        } // done adjusting the damping factor

      } // end of iterative minimization

      // if reached max # of iterations throw if the error is too terrible
      assert(not (iter == maxniter && errnormI > 1e-10));

      for(unsigned int i=0; i<n; ++i)
        geminal[i] = std::make_pair(aa[i], cc[i]);
    }
#endif

} // end of namespace libint2

#endif // C++ only
#endif // header guard
