/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
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

#ifndef _libint2_src_lib_libint_boys_h_
#define _libint2_src_lib_libint_boys_h_

#if defined(__cplusplus)

#include <libint2/util/optional.h>
#include <libint2/util/vector.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef _MSC_VER
#define posix_memalign(p, a, s) \
  (((*(p)) = _aligned_malloc((s), (a))), *(p) ? 0 : errno)
#define posix_memfree(p) ((_aligned_free((p))))
#else
#define posix_memfree(p) ((free((p))))
#endif

// from now on at least C++11 is required by default
#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "Libint2 C++ API requires C++11 support"
#endif

#include <libint2/boys_fwd.h>
#include <libint2/initialize.h>
#include <libint2/numeric.h>

#if HAVE_LAPACK  // use F77-type interface for now, switch to LAPACKE later
extern "C" void dgesv_(const int* n, const int* nrhs, double* A, const int* lda,
                       int* ipiv, double* b, const int* ldb, int* info);
#endif

namespace libint2 {

/// holds tables of expensive quantities
template <typename Real>
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
      if (idf >= 1) df[1] = 1.0;
      if (idf >= 2) df[2] = 1.0;
      for (int i = 3; i <= idf; i++) {
        df[i] = (i - 1) * df[i - 2];
      }
    }

    if (ibc >= 0) {
      bc_.resize((ibc + 1) * (ibc + 1));
      std::fill(bc_.begin(), bc_.end(), Real(0));
      bc.resize(ibc + 1);
      bc[0] = &bc_[0];
      for (int i = 1; i <= ibc; ++i) bc[i] = bc[i - 1] + (ibc + 1);

      for (int i = 0; i <= ibc; i++) bc[i][0] = 1.0;
      for (int i = 0; i <= ibc; i++)
        for (int j = 1; j <= i; ++j)
          bc[i][j] = bc[i][j - 1] * Real(i - j + 1) / Real(j);
    }

    for (int i = 0; i < 128; i++) {
      twoi1[i] = 1.0 / (Real(2.0) * i + Real(1.0));
      ihalf[i] = Real(i) - Real(0.5);
    }
  }

  ~ExpensiveNumbers() {}

  std::vector<Real> fac;
  std::vector<Real> df;
  std::vector<Real*> bc;

  // these quantitites are needed with indices <= mmax
  // 64 is sufficient to handle up to 4 center integrals with up to L=15 basis
  // functions but need higher values for Yukawa integrals ...
  Real twoi1[128]; /* 1/(2 i + 1); needed for downward recursion */
  Real ihalf[128]; /* i - 0.5, needed for upward recursion */

 private:
  std::vector<Real> bc_;
};

/** Computes the Boys function, \f$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \,
 * {\rm d}u \f$, using single algorithm (asymptotic expansion). Slow for the
 * sake of precision control. Useful in two cases: <ul> <li> for reference
 * purposes, if \c Real supports high/arbitrary precision, and </li> <li> for
 * moderate values of \f$ T \f$, if \c Real is a low-precision floating-point
 * type. N.B. FmEval_Reference2 , which can compute for all practical values of
 * \f$ T \f$ and \f$ m \f$, is recommended with standard \c Real types (\c
 * double and \c float). </li>
 * </ul>
 *
 * \note Precision is controlled heuristically, i.e. cannot be guaranteed
 * mathematically; will stop if absolute precision is reached, or precision of
 * \c Real is exhausted. It is important that \c std::numeric_limits<Real> is
 * defined appropriately.
 *
 * @tparam Real the type to use for all floating-point computations.
 *         Must be able to compute logarithm and exponential, i.e.
 *         log(x) and exp(x), where x is Real, must be valid expressions.
 */
template <typename Real>
struct FmEval_Reference {
  static std::shared_ptr<const FmEval_Reference> instance(
      int /* mmax */, Real /* precision */) {
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ = std::make_shared<const FmEval_Reference>();

    return instance_;
  }

  /// computes a single value of \f$ F_m(T) \f$ using MacLaurin series to full
  /// precision of @c Real
  static Real eval(Real T, size_t m) {
    assert(m < 100);
    static const Real T_crit =
        std::numeric_limits<Real>::is_bounded
            ? -log(std::numeric_limits<Real>::min() * 100.5 / 2.)
            : Real(0);
    if (std::numeric_limits<Real>::is_bounded && T > T_crit)
      throw std::overflow_error(
          "FmEval_Reference<Real>::eval: Real lacks precision for the given "
          "value of argument T");
    static const Real half = Real(1) / 2;
    Real denom = (m + half);
    using std::exp;
    Real term = exp(-T) / (2 * denom);
    Real old_term = 0;
    Real sum = term;
    const Real epsilon = get_epsilon(T);
    const Real epsilon_divided_10 = epsilon / 10;
    do {
      denom += 1;
      old_term = term;
      term = old_term * T / denom;
      sum += term;
      // rel_error = term / sum , hence iterate until rel_error = epsilon
      //  however, must ensure that contributions are decreasing to ensure that
      //  omitted contributions are smaller than epsilon
    } while (term > sum * epsilon_divided_10 || old_term < term);

    return sum;
  }

  /// fills up an array of Fm(T) for m in [0,mmax]
  /// @param[out] Fm array to be filled in with the Boys function values, must
  /// be at least mmax+1 elements long
  /// @param[in] T the Boys function argument
  /// @param[in] mmax the maximum value of m for which Boys function will be
  /// computed;
  static void eval(Real* Fm, Real T, size_t mmax) {
    // evaluate for mmax using MacLaurin series
    // it converges fastest for the largest m -> use it to compute Fmmax(T)
    //  see JPC 94, 5564 (1990).
    for (size_t m = 0; m <= mmax; ++m) Fm[m] = eval(T, m);
    return;
  }
};

/** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \,
 * {\rm d}u \$, using multi-algorithm approach (upward recursion for T>=117, and
 * asymptotic summation for T<117). This is slow and should be used for
 * reference purposes, e.g. computing the interpolation tables. Precision is not
 * always guaranteed as it is limited by the precision of \c Real type. When \c
 * Real is \c double, can maintain absolute precision of epsilon for up to m=40.
 *
 * @tparam Real the type to use for all floating-point computations.
 *         Must be able to compute logarithm, exponential, square root, and
 * error function, i.e. log(x), exp(x), sqrt(x), and erf(x), where x is Real,
 * must be valid expressions.
 */
template <typename Real>
struct FmEval_Reference2 {
  static std::shared_ptr<const FmEval_Reference2> instance(
      int /* mmax */, Real /* precision */) {
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ = std::make_shared<const FmEval_Reference2>();

    return instance_;
  }

  /// fills up an array of Fm(T) for m in [0,mmax]
  /// @param[out] Fm array to be filled in with the Boys function values, must
  /// be at least mmax+1 elements long
  /// @param[in] t the Boys function argument
  /// @param[in] mmax the maximum value of m for which Boys function will be
  /// computed;
  static void eval(Real* Fm, Real t, size_t mmax) {
    if (t < Real(117)) {
      FmEval_Reference<Real>::eval(Fm, t, mmax);
    } else {
      const Real two_over_sqrt_pi{
          1.128379167095512573896158903121545171688101258657997713688171443421284936882986828973487320404214727};
      const Real K = 1 / two_over_sqrt_pi;

      auto t2 = 2 * t;
      auto et = exp(-t);
      auto sqrt_t = sqrt(t);
      Fm[0] = K * erf(sqrt_t) / sqrt_t;
      if (mmax > 0)
        for (size_t m = 0; m <= mmax - 1; m++) {
          Fm[m + 1] = ((2 * m + 1) * Fm[m] - et) / (t2);
        }
    }
  }
};

/** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \,
 * {\rm d}u \$, using 7-th order Chebyshev interpolation.
 */
template <typename Real = double>
class FmEval_Chebyshev7 {
#include <libint2/boys_cheb7.h>

  static_assert(std::is_same<Real, double>::value,
                "FmEval_Chebyshev7 only supports double as the real type");

  static constexpr int ORDER = interpolation_order;  //!, interpolation order
  static constexpr int ORDERp1 = ORDER + 1;          //!< ORDER + 1

  static constexpr Real T_crit =
      cheb_table_tmax;  //!< critical value of T above which safe to use upward
                        //!< recursion
  static constexpr Real delta = cheb_table_delta;    //!< interval size
  static constexpr Real one_over_delta = 1 / delta;  //! 1/delta

  int mmax;  //!< the maximum m that is tabulated
  ExpensiveNumbers<double> numbers_;
  Real* c; /* the Chebyshev coefficients table, T_crit*one_over_delta by
              mmax*ORDERp1 */

 public:
  /// \param m_max maximum value of the Boys function index; set to -1 to skip
  /// initialization \param precision the desired relative precision
  /// \throw std::invalid_argument if \c m_max is greater than
  ///   \c cheb_table_mmax (see boys_cheb7.h)
  /// \throw std::invalid_argument if \c precision is smaller than
  ///   `std::numeric_limits<double>::epsilon()`
  FmEval_Chebyshev7(int m_max,
                    double precision = std::numeric_limits<double>::epsilon())
      : mmax(m_max), numbers_(14) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || \
    defined(_M_IX86)
#if !defined(__AVX__) && defined(NDEBUG)
    if (libint2::verbose()) {
      static bool printed_performance_warning = false;
      if (!printed_performance_warning) {
        libint2::verbose_stream()
            << "libint2::FmEval_Chebyshev7 on x86(-64) platforms needs AVX "
               "support for best performance"
            << std::endl;
        printed_performance_warning = true;
      }
    }
#endif
#endif
    if (precision < std::numeric_limits<double>::epsilon())
      throw std::invalid_argument(
          std::string(
              "FmEval_Chebyshev7 does not support precision smaller than ") +
          std::to_string(std::numeric_limits<double>::epsilon()));
    if (mmax > cheb_table_mmax)
      throw std::invalid_argument(
          "FmEval_Chebyshev7::init() : requested mmax exceeds the "
          "hard-coded mmax");
    if (m_max >= 0) init_table();
  }
  ~FmEval_Chebyshev7() {
    if (mmax >= 0) {
      posix_memfree(c);
    }
  }

  /// Singleton interface allows to manage the lone instance; adjusts max m
  /// values as needed in thread-safe fashion
  static std::shared_ptr<const FmEval_Chebyshev7> instance(int m_max,
                                                           double = 0.0) {
    assert(m_max >= 0);
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ = std::make_shared<const FmEval_Chebyshev7>(m_max);

    while (instance_->max_m() < m_max) {
      static std::mutex mtx;
      std::lock_guard<std::mutex> lck(mtx);
      if (instance_->max_m() < m_max) {
        auto new_instance = std::make_shared<const FmEval_Chebyshev7>(m_max);
        instance_ = new_instance;
      }
    }

    return instance_;
  }

  /// @return the maximum value of m for which the Boys function can be computed
  /// with this object
  int max_m() const { return mmax; }

  /// fills in Fm with computed Boys function values for m in [0,mmax]
  /// @param[out] Fm array to be filled in with the Boys function values, must
  /// be at least mmax+1 elements long
  /// @param[in] x the Boys function argument
  /// @param[in] mmax the maximum value of m for which Boys function will be
  /// computed; mmax must be <= the value returned by max_m
  inline void eval(Real* Fm, Real x, int m_max) const {
    // large T => use upward recursion
    // cost = 1 div + 1 sqrt + (1 + 2*(m-1)) muls
    if (x > T_crit) {
      const double one_over_x = 1 / x;
      Fm[0] = 0.88622692545275801365 *
              sqrt(one_over_x);  // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
      if (m_max == 0) return;
      // this upward recursion formula omits - e^(-x)/(2x), which for x>T_crit
      // is small enough to guarantee full double precision
      for (int i = 1; i <= m_max; i++)
        Fm[i] = Fm[i - 1] * numbers_.ihalf[i] * one_over_x;  // see Eq. (9.8.13)
      return;
    }

    // ---------------------------------------------
    // small and intermediate arguments => interpolate Fm and (optional)
    // downward recursion
    // ---------------------------------------------
    // which interval does this x fall into?
    const Real x_over_delta = x * one_over_delta;
    const int iv = int(x_over_delta);  // the interval index
    const Real xd =
        x_over_delta - (Real)iv - 0.5;  // this ranges from -0.5 to 0.5
    const int m_min = 0;

#if defined(__AVX__)
    const auto x2 = xd * xd;
    const auto x3 = x2 * xd;
    const auto x4 = x2 * x2;
    const auto x5 = x2 * x3;
    const auto x6 = x3 * x3;
    const auto x7 = x3 * x4;
    libint2::simd::VectorAVXDouble x0vec(1., xd, x2, x3);
    libint2::simd::VectorAVXDouble x1vec(x4, x5, x6, x7);
#endif  // AVX

    const Real* d =
        c + (ORDERp1) * (iv * (mmax + 1) +
                         m_min);  // ptr to the interpolation data for m=mmin
    int m = m_min;
#if defined(__AVX__)
    if (m_max - m >= 3) {
      const int unroll_size = 4;
      const int m_fence = (m_max + 2 - unroll_size);
      for (; m < m_fence; m += unroll_size, d += ORDERp1 * unroll_size) {
        libint2::simd::VectorAVXDouble d00v, d01v, d10v, d11v, d20v, d21v, d30v,
            d31v;
        d00v.load_aligned(d);
        d01v.load_aligned((d + 4));
        d10v.load_aligned(d + ORDERp1);
        d11v.load_aligned((d + 4) + ORDERp1);
        d20v.load_aligned(d + 2 * ORDERp1);
        d21v.load_aligned((d + 4) + 2 * ORDERp1);
        d30v.load_aligned(d + 3 * ORDERp1);
        d31v.load_aligned((d + 4) + 3 * ORDERp1);
        libint2::simd::VectorAVXDouble fm0 = d00v * x0vec + d01v * x1vec;
        libint2::simd::VectorAVXDouble fm1 = d10v * x0vec + d11v * x1vec;
        libint2::simd::VectorAVXDouble fm2 = d20v * x0vec + d21v * x1vec;
        libint2::simd::VectorAVXDouble fm3 = d30v * x0vec + d31v * x1vec;
        libint2::simd::VectorAVXDouble sum0123 =
            horizontal_add(fm0, fm1, fm2, fm3);
        sum0123.convert(&Fm[m]);
      }
    }  // unroll_size=4
    if (m_max - m >= 1) {
      const int unroll_size = 2;
      const int m_fence = (m_max + 2 - unroll_size);
      for (; m < m_fence; m += unroll_size, d += ORDERp1 * unroll_size) {
        libint2::simd::VectorAVXDouble d00v, d01v, d10v, d11v;
        d00v.load_aligned(d);
        d01v.load_aligned((d + 4));
        d10v.load_aligned(d + ORDERp1);
        d11v.load_aligned((d + 4) + ORDERp1);
        libint2::simd::VectorAVXDouble fm0 = d00v * x0vec + d01v * x1vec;
        libint2::simd::VectorAVXDouble fm1 = d10v * x0vec + d11v * x1vec;
        libint2::simd::VectorSSEDouble sum01 = horizontal_add(fm0, fm1);
        sum01.convert(&Fm[m]);
      }
    }  // unroll_size=2
    {  // no unrolling
      for (; m <= m_max; ++m, d += ORDERp1) {
        libint2::simd::VectorAVXDouble d0v, d1v;
        d0v.load_aligned(d);
        d1v.load_aligned(d + 4);
        Fm[m] = horizontal_add(d0v * x0vec + d1v * x1vec);
      }
    }
#else  // AVX not available
    for (m = m_min; m <= m_max; ++m, d += ORDERp1) {
      Fm[m] =
          d[0] +
          xd * (d[1] +
                xd * (d[2] +
                      xd * (d[3] +
                            xd * (d[4] +
                                  xd * (d[5] + xd * (d[6] + xd * (d[7])))))));

      //        // check against the reference value
      //        if (false) {
      //          double refvalue = FmEval_Reference2<double>::eval(x, m); //
      //          compute F(T) if (abs(refvalue - Fm[m]) > 1e-10) {
      //            std::cout << "T = " << x << " m = " << m << " cheb = "
      //                      << Fm[m] << " ref = " << refvalue << std::endl;
      //          }
      //        }
    }
#endif

  }  // eval()

 private:
  void init_table() {
    // get memory
    void* result;
    int status = posix_memalign(
        &result, ORDERp1 * sizeof(Real),
        (mmax + 1) * cheb_table_nintervals * ORDERp1 * sizeof(Real));
    if (status != 0) {
      if (status == EINVAL)
        throw std::logic_error(
            "FmEval_Chebyshev7::init() : posix_memalign failed, alignment must "
            "be a power of 2 at least as large as sizeof(void *)");
      if (status == ENOMEM) throw std::bad_alloc();
      abort();  // should be unreachable
    }
    c = static_cast<Real*>(result);

    // copy contents of static table into c
    // need all intervals
    for (std::size_t iv = 0; iv < cheb_table_nintervals; ++iv) {
      // but only values of m up to mmax
      std::copy(cheb_table[iv], cheb_table[iv] + (mmax + 1) * ORDERp1,
                c + (iv * (mmax + 1)) * ORDERp1);
    }
  }

};  // FmEval_Chebyshev7

#if LIBINT2_CONSTEXPR_STATICS
template <typename Real>
constexpr double FmEval_Chebyshev7<
    Real>::cheb_table[FmEval_Chebyshev7<Real>::cheb_table_nintervals]
                     [(FmEval_Chebyshev7<Real>::cheb_table_mmax + 1) *
                      (FmEval_Chebyshev7<Real>::interpolation_order + 1)];
#else
// clang needs an explicit specifalization declaration to avoid warning
#ifdef __clang__
template <typename Real>
double FmEval_Chebyshev7<
    Real>::cheb_table[FmEval_Chebyshev7<Real>::cheb_table_nintervals]
                     [(FmEval_Chebyshev7<Real>::cheb_table_mmax + 1) *
                      (FmEval_Chebyshev7<Real>::interpolation_order + 1)];
#endif
#endif

#ifndef STATIC_OON
#define STATIC_OON
namespace {
const double oon[] = {0.0,       1.0,       1.0 / 2.0,  1.0 / 3.0,
                      1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0,  1.0 / 7.0,
                      1.0 / 8.0, 1.0 / 9.0, 1.0 / 10.0, 1.0 / 11.0};
}
#endif

/** Computes the Boys function, \$ F_m (T) = \int_0^1 u^{2m} \exp(-T u^2) \,
 * {\rm d}u \$, using Taylor interpolation of up to 8-th order.
 * @tparam Real the type to use for all floating-point computations. Following
 * expressions must be valid: exp(Real), pow(Real), fabs(Real), max(Real), and
 * floor(Real).
 * @tparam INTERPOLATION_ORDER the interpolation order. The higher the order the
 * less memory this object will need, but the computational cost will increase
 * (usually very slightly)
 */
template <typename Real = double, int INTERPOLATION_ORDER = 7>
class FmEval_Taylor {
 public:
  static_assert(std::is_same<Real, double>::value,
                "FmEval_Taylor only supports double as the real type");

  static const int max_interp_order = 8;
  static const bool INTERPOLATION_AND_RECURSION =
      false;  // compute F_lmax(T) and then iterate down to F_0(T)? Else use
              // interpolation only
  const double soft_zero_;

  /// Constructs the object to be able to compute Boys funcion for m in
  /// [0,mmax], with relative \c precision
  FmEval_Taylor(unsigned int mmax,
                Real precision = std::numeric_limits<Real>::epsilon())
      : soft_zero_(1e-6),
        cutoff_(precision),
        numbers_(INTERPOLATION_ORDER + 1, 2 * (mmax + INTERPOLATION_ORDER)) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || \
    defined(_M_IX86)
#if !defined(__AVX__) && defined(NDEBUG)
    if (libint2::verbose()) {
      static bool printed_performance_warning = false;
      if (!printed_performance_warning) {
        libint2::verbose_stream()
            << "libint2::FmEval_Taylor on x86(-64) platforms needs AVX support "
               "for best performance"
            << std::endl;
        printed_performance_warning = true;
      }
    }
#endif
#endif

    assert(mmax <= 63);

    const double sqrt_pi = std::sqrt(M_PI);

    /*---------------------------------------
     We are doing Taylor interpolation with
     n=TAYLOR_ORDER terms here:
     error <= delT^n/(n+1)!
     ---------------------------------------*/
    delT_ = 2.0 * std::pow(cutoff_ * numbers_.fac[INTERPOLATION_ORDER + 1],
                           1.0 / INTERPOLATION_ORDER);
    oodelT_ = 1.0 / delT_;
    max_m_ = mmax + INTERPOLATION_ORDER;

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
      const double egamma =
          cutoff_ * sqrt_pi * numbers_.df[2 * m] / std::pow(2.0, m);
      double T_new = T;
      double func;
      do {
        const double damping_factor = 0.2;
        T = T_new;
        /* f(T) = the difference between LHS and RHS of the equation above */
        func = std::pow(T, m - 0.5) * std::exp(-T) - egamma;
        const double dfuncdT =
            ((m - 0.5) * std::pow(T, m - 1.5) - std::pow(T, m - 0.5)) *
            std::exp(-T);
        /* f(T) has 2 roots and has a maximum in between. If f'(T) > 0 we are to
         * the left of the hump. Make a big step to the right. */
        if (dfuncdT > 0.0) {
          T_new *= 2.0;
        } else {
          /* damp the step */
          double deltaT = -func / dfuncdT;
          const double sign_deltaT = (deltaT > 0.0) ? 1.0 : -1.0;
          const double max_deltaT = damping_factor * T;
          if (std::fabs(deltaT) > max_deltaT) deltaT = sign_deltaT * max_deltaT;
          T_new = T + deltaT;
        }
        if (T_new <= 0.0) {
          T_new = T / 2.0;
        }
      } while (std::fabs(func / egamma) >= soft_zero_);
      T_crit_[m] = T_new;
      const int T_idx = (int)std::floor(T_new / delT_);
      max_T_ = std::max(max_T_, T_idx);
    }

    // allocate the grid (see the comments below)
    {
      const int nrow = max_T_ + 1;
      const int ncol = max_m_ + 1;
      grid_ = new Real*[nrow];
      grid_[0] = new Real[nrow * ncol];
      // std::cout << "Allocated interpolation table of " << nrow * ncol << "
      // reals" << std::endl;
      for (int r = 1; r < nrow; ++r) grid_[r] = grid_[r - 1] + ncol;
    }

    /*-------------------------------------------------------
     Tabulate the gamma function from t=delT to T_crit[m]:
     1) include T=0 though the table is empty for T=0 since
     Fm(0) is simple to compute
     -------------------------------------------------------*/
    /*--- do the mmax first ---*/
    for (int T_idx = max_T_; T_idx >= 0; --T_idx) {
      const double T = T_idx * delT_;
      libint2::FmEval_Reference2<double>::eval(grid_[T_idx], T, max_m_);
    }
  }

  ~FmEval_Taylor() {
    delete[] T_crit_;
    delete[] grid_[0];
    delete[] grid_;
  }

  /// Singleton interface allows to manage the lone instance;
  /// adjusts max m and precision values as needed in thread-safe fashion
  static std::shared_ptr<const FmEval_Taylor> instance(
      unsigned int mmax,
      Real precision = std::numeric_limits<Real>::epsilon()) {
    assert(mmax >= 0);
    assert(precision >= 0);
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ =
        std::make_shared<const FmEval_Taylor>(mmax, precision);

    while (instance_->max_m() < mmax || instance_->precision() > precision) {
      static std::mutex mtx;
      std::lock_guard<std::mutex> lck(mtx);
      if (instance_->max_m() < mmax || instance_->precision() > precision) {
        auto new_instance =
            std::make_shared<const FmEval_Taylor>(mmax, precision);
        instance_ = new_instance;
      }
    }

    return instance_;
  }

  /// @return the maximum value of m for which this object can compute the Boys
  /// function
  int max_m() const { return max_m_ - INTERPOLATION_ORDER + 1; }
  /// @return the precision with which this object can compute the Boys function
  Real precision() const { return cutoff_; }

  /// computes Boys function values with m index in range [0,mmax]
  /// @param[out] Fm array to be filled in with the Boys function values, must
  /// be at least mmax+1 elements long
  /// @param[in] x the Boys function argument
  /// @param[in] mmax the maximum value of m for which Boys function will be
  /// computed;
  ///                  it must be <= the value returned by max_m() (this is not
  ///                  checked)
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
        const double one_over_x = 1.0 / T;
        Fm[0] =
            0.88622692545275801365 *
            sqrt(one_over_x);  // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
        if (mmax == 0) return;
        // this upward recursion formula omits - e^(-x)/(2x), which for x>T_crit
        // is <1e-15
        for (int i = 1; i <= mmax; i++)
          Fm[i] =
              Fm[i - 1] * numbers_.ihalf[i] * one_over_x;  // see Eq. (9.8.13)
        return;
      }
    }

    // since Tcrit grows with mmax, this condition only needs to be determined
    // once
    if (T > T_crit_[mmax]) {
      double pow_two_T_to_minusjp05 = std::pow(two_T, -mmax - 0.5);
      for (int m = mmax; m >= mmin; --m) {
        /*--- Asymptotic formula ---*/
        Fm[m] = numbers_.df[2 * m] * sqrt_pio2 * pow_two_T_to_minusjp05;
        pow_two_T_to_minusjp05 *= two_T;
      }
    } else {
      const int T_ind = (int)(0.5 + T * oodelT_);
      const Real h = T_ind * delT_ - T;
      const Real* F_row = grid_[T_ind] + mmin;

#if defined(__AVX__)
      libint2::simd::VectorAVXDouble h0123, h4567;
      if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
        const double h2 = h * h * oon[2];
        const double h3 = h2 * h * oon[3];
        h0123 = libint2::simd::VectorAVXDouble(1.0, h, h2, h3);
        if (INTERPOLATION_ORDER == 7) {
          const double h4 = h3 * h * oon[4];
          const double h5 = h4 * h * oon[5];
          const double h6 = h5 * h * oon[6];
          const double h7 = h6 * h * oon[7];
          h4567 = libint2::simd::VectorAVXDouble(h4, h5, h6, h7);
        }
      }
      //          libint2::simd::VectorAVXDouble h0123(1.0);
      //          libint2::simd::VectorAVXDouble h4567(1.0);
#endif

      int m = mmin;
      if (mmax - m >= 1) {
        const int unroll_size = 2;
        const int m_fence = (mmax + 2 - unroll_size);
        for (; m < m_fence; m += unroll_size, F_row += unroll_size) {
#if defined(__AVX__)
          if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
            libint2::simd::VectorAVXDouble fr0_0123;
            fr0_0123.load(F_row);
            libint2::simd::VectorAVXDouble fr1_0123;
            fr1_0123.load(F_row + 1);
            libint2::simd::VectorSSEDouble fm01 =
                horizontal_add(fr0_0123 * h0123, fr1_0123 * h0123);
            if (INTERPOLATION_ORDER == 7) {
              libint2::simd::VectorAVXDouble fr0_4567;
              fr0_4567.load(F_row + 4);
              libint2::simd::VectorAVXDouble fr1_4567;
              fr1_4567.load(F_row + 5);
              fm01 += horizontal_add(fr0_4567 * h4567, fr1_4567 * h4567);
            }
            fm01.convert(&Fm[m]);
          } else {
#endif
            Real total0 = 0.0, total1 = 0.0;
            for (int i = INTERPOLATION_ORDER; i >= 1; --i) {
              total0 = oon[i] * h * (F_row[i] + total0);
              total1 = oon[i] * h * (F_row[i + 1] + total1);
            }
            Fm[m] = F_row[0] + total0;
            Fm[m + 1] = F_row[1] + total1;
#if defined(__AVX__)
          }
#endif
        }
      }                 // unroll_size = 2
      if (m <= mmax) {  // unroll_size = 1
#if defined(__AVX__)
        if (INTERPOLATION_ORDER == 3 || INTERPOLATION_ORDER == 7) {
          libint2::simd::VectorAVXDouble fr0123;
          fr0123.load(F_row);
          if (INTERPOLATION_ORDER == 7) {
            libint2::simd::VectorAVXDouble fr4567;
            fr4567.load(F_row + 4);
            //                libint2::simd::VectorSSEDouble fm =
            //                horizontal_add(fr0123*h0123, fr4567*h4567); Fm[m]
            //                = horizontal_add(fm);
            Fm[m] = horizontal_add(fr0123 * h0123 + fr4567 * h4567);
          } else {  // INTERPOLATION_ORDER == 3
            Fm[m] = horizontal_add(fr0123 * h0123);
          }
        } else {
#endif
          Real total = 0.0;
          for (int i = INTERPOLATION_ORDER; i >= 1; --i) {
            total = oon[i] * h * (F_row[i] + total);
          }
          Fm[m] = F_row[0] + total;
#if defined(__AVX__)
        }
#endif
      }  // unroll_size = 1

      // check against the reference value
      //          if (false) {
      //            double refvalue = FmEval_Reference2<double>::eval(T, mmax);
      //            // compute F(T) with m=mmax if (abs(refvalue - Fm[mmax]) >
      //            1e-14) {
      //              std::cout << "T = " << T << " m = " << mmax << " cheb = "
      //                  << Fm[mmax] << " ref = " << refvalue << std::endl;
      //            }
      //          }

    }  // if T < T_crit

    /*------------------------------------
     And then do downward recursion in j
     ------------------------------------*/
    if (INTERPOLATION_AND_RECURSION && mmin > 0) {
      const Real exp_mT = std::exp(-T);
      for (int m = mmin - 1; m >= 0; --m) {
        Fm[m] = (exp_mT + two_T * Fm[m + 1]) * numbers_.twoi1[m];
      }
    }
  }

 private:
  Real** grid_;  /* Table of "exact" Fm(T) values. Row index corresponds to
    values of T (max_T+1 rows), column index to values
    of m (max_m+1 columns) */
  Real delT_;    /* The step size for T, depends on cutoff */
  Real oodelT_;  /* 1.0 / delT_, see above */
  Real cutoff_;  /* Tolerance cutoff used in all computations of Fm(T) */
  int max_m_;    /* Maximum value of m in the table, depends on cutoff
      and the number of terms in Taylor interpolation */
  int max_T_;    /* Maximum index of T in the table, depends on cutoff
      and m */
  Real* T_crit_; /* Maximum T for each row, depends on cutoff;
   for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
   for a given m and T_idx > max_T_idx[m] use the asymptotic formula */

  ExpensiveNumbers<double> numbers_;

  /**
   * Power series estimate of the error introduced by replacing
   * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with
   * analytically integrable \f$ \int_0^\infty \exp(-T t^2) t^{2 m} \,
   * \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
   * @param m
   * @param T
   * @return the error estimate
   */
  static double truncation_error(unsigned int m, double T) {
    const double m2 = m * m;
    const double m3 = m2 * m;
    const double m4 = m2 * m2;
    const double T2 = T * T;
    const double T3 = T2 * T;
    const double T4 = T2 * T2;
    const double T5 = T2 * T3;

    const double result =
        exp(-T) *
        (105 + 16 * m4 + 16 * m3 * (T - 8) - 30 * T + 12 * T2 - 8 * T3 +
         16 * T4 + 8 * m2 * (43 - 9 * T + 2 * T2) +
         4 * m * (-88 + 23 * T - 8 * T2 + 4 * T3)) /
        (32 * T5);
    return result;
  }
  /**
   * Leading-order estimate of the error introduced by replacing
   * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with
   * analytically integrable \f$ \int_0^\infty \exp(-T t^2) t^{2 m} \,
   * \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
   * @param m
   * @param T
   * @return the error estimate
   */
  static double truncation_error(double T) {
    const double result = exp(-T) / (2 * T);
    return result;
  }
};

namespace detail {
constexpr double pow10(long exp) {
  return (exp == 0) ? 1.
                    : ((exp > 0) ? 10. * pow10(exp - 1) : 0.1 * pow10(exp + 1));
}
}  // namespace detail

//////////////////////////////////////////////////////////
/// core integral for Yukawa and exponential interactions
//////////////////////////////////////////////////////////

/**
 * Evaluates core integral for Gaussian integrals over the Yukawa potential
 * \f$ \exp(- \zeta r) / r \f$ and the exponential interaction
 * \f$ \exp(- \zeta r) \f$
 * @tparam Real real type
 * @warning only @p Real = double is supported
 * @warning guarantees absolute precision of only about 1e-14
 * @warning for vast majority of practical uses the integral is computed
 * by interpolation or closed-form evaluation. But rarely the
 * the "analytic" route will cause overflow, so will represent kernel
 * as a linear combination of Gaussians to compute without overflow.
 */
template <typename Real = double>
struct TennoGmEval {
 private:
#include <libint2/tenno_cheb.h>

  static_assert(std::is_same<Real, double>::value,
                "TennoGmEval only supports double as the real type");

  static const int mmin_ = -1;
  static constexpr Real Tmax =
      (1 << cheb_table_tmaxlog2);  //!< critical value of T above which use
                                   //!< upward recursion
  static constexpr Real Umax = detail::pow10(
      cheb_table_umaxlog10);  //!< max value of U for which to interpolate
  static constexpr Real Umin = detail::pow10(
      cheb_table_uminlog10);  //!< min value of U for which to interpolate
  static constexpr Real one_over_Umin =
      1. / Umin;  //!< min value of U for which to interpolate
  static constexpr std::size_t ORDERp1 = interpolation_order + 1;
  static constexpr Real maxabsprecision =
      1.4e-14;  //!< guaranteed abs precision of the interpolation table for m>0
                // map of zeta to existing gaussian fit to slater to avoid
                // recomputing

  // to avoid recomputing Gaussian fits of the exponential cache
  // them for every value of zeta
  using gaussian_fit_t = std::vector<std::pair<double, double>>;
  using gaussian_fits_t = std::map<double, gaussian_fit_t>;  // zeta -> fit
  // TODO replace by atomic shared_ptr to avoid leaking memory at the end
  using gaussian_fits_aptr_t = std::atomic<gaussian_fits_t*>;
  static gaussian_fits_aptr_t& gaussian_fits_aptr_accessor() {
    static gaussian_fits_aptr_t gaussian_fits_aptr_;
    return gaussian_fits_aptr_;
  }
  static gaussian_fit_t& get_gaussian_fit(double zeta) {
    while (gaussian_fits_aptr_accessor().load() == nullptr ||
           gaussian_fits_aptr_accessor().load()->count(zeta) == 0) {
      // remember ptr to the cache we will be adding new fit to
      auto* gaussian_fits_ptr = gaussian_fits_aptr_accessor().load();

      // use existing cache, or create a new one
      auto new_gaussian_fits_ptr = new gaussian_fits_t(
          gaussian_fits_ptr == nullptr ? gaussian_fits_t()
                                       : *gaussian_fits_ptr);
      // constexpr double low_end_range = 1e-5;
      // constexpr double high_end_range = 5;
      // constexpr double eps = 1e-12;
      std::vector<std::pair<double, double>> geminal =
          fit_slater(zeta /* , low_end_range, high_end_range, eps */);
      gaussian_fits_t::iterator it;
      bool inserted;
      std::tie(it, inserted) = new_gaussian_fits_ptr->emplace(zeta, geminal);

      // inserted may be false if someone inserted between count query and
      // grabbing gaussian_fits_ptr ... if inserted new fit try to replace the
      // value in gaussian_fits_aptr_
      if (inserted) {
        auto* expected_gaussian_fits_ptr = gaussian_fits_ptr;
        if (gaussian_fits_aptr_accessor().compare_exchange_strong(
                expected_gaussian_fits_ptr, new_gaussian_fits_ptr)) {
          // can't delete since someone might be using it ... need atomic
          // shared_ptr
          // delete gaussian_fits_ptr;
          break;
        } else {
          delete new_gaussian_fits_ptr;
        }
      }
    }
    return gaussian_fits_aptr_accessor().load()->at(zeta);
  }
  mutable std::unique_ptr<GaussianGmEval<Real, -1>> yukawa_gaussian_fit_eval_;
  mutable std::unique_ptr<GaussianGmEval<Real, 0>> slater_gaussian_fit_eval_;

 public:
  /// \param m_max maximum value of the Gm function index
  /// \param precision the desired *absolute* precision (relative precision for
  /// most intervals will be below epsilon, but for large T/U values and high m
  /// relative precision is low); since evaluation is done primarily
  /// by interpolation, precision is not controlled dynamically, thus this
  /// should not be used.
  /// \throw std::invalid_argument if \c m_max is greater than
  ///   \c cheb_table_mmax (see tenno_cheb.h)
  /// \throw std::invalid_argument if \c precision is smaller
  ///   than \c maxabsprecision
  TennoGmEval(unsigned int mmax, Real precision = maxabsprecision)
      : mmax_(mmax), precision_(precision), numbers_() {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || \
    defined(_M_IX86)
#if !defined(__AVX__) && defined(NDEBUG)
    if (libint2::verbose()) {
      static bool printed_performance_warning = false;
      if (!printed_performance_warning) {
        libint2::verbose_stream()
            << "libint2::TennoGmEval on x86(-64) platforms needs AVX support "
               "for best performance"
            << std::endl;
        printed_performance_warning = true;
      }
    }
#endif
#endif

    if (precision_ < maxabsprecision)
      throw std::invalid_argument(
          std::string("TennoGmEval does not support precision smaller than ") +
          std::to_string(maxabsprecision));

    if (mmax > cheb_table_mmax)
      throw std::invalid_argument(
          "TennoGmEval::init() : requested mmax exceeds the "
          "hard-coded mmax");
    init_table();
  }

  ~TennoGmEval() {
    if (c_ != nullptr) posix_memfree(c_);
  }

  /// Singleton interface allows to manage the lone instance; adjusts max m
  /// values as needed in thread-safe fashion
  static std::shared_ptr<const TennoGmEval> instance(int m_max, double = 0) {
    assert(m_max >= 0);
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ = std::make_shared<const TennoGmEval>(m_max);

    while (instance_->max_m() < m_max) {
      static std::mutex mtx;
      std::lock_guard<std::mutex> lck(mtx);
      if (instance_->max_m() < m_max) {
        auto new_instance = std::make_shared<const TennoGmEval>(m_max);
        instance_ = new_instance;
      }
    }

    return instance_;
  }

  unsigned int max_m() const { return mmax_; }
  /// @return the precision with which this object can compute the result
  Real precision() const { return precision_; }

  /// @param[in] Gm pointer to array of @c mmax+1 @c Real elements, on
  ///               return this contains the core integral for Yukawa
  ///               interaction, \f$ \exp(-zeta r_{12})/r_{12} \f$ ,
  ///               namely \f$ G_{m}(T,U) = \int_0^1 t^{2m} \exp(U(1-t^{-2}) -
  ///               Tt^2) dt, m \in [0,m_\mathrm{max}] \f$, where \f$ T = \rho
  ///               |\vec{P} - \vec{Q}|^2 \f$ and \f$ U = \zeta^2 / (4 \rho) \f$
  /// @param[in] one_over_rho \f$ 1/\rho \f$
  /// @param[in] T \f$ T \f$
  /// @param[in] mmax \f$ m_\mathrm{max} \f$
  /// @param[in] zeta \f$ \zeta \f$
  /// @throw std::invalid_argument if \f$ \zeta^2/(4 \rho) \f$ exceeds \c Umax
  void eval_yukawa(Real* Gm, Real one_over_rho, Real T, size_t mmax,
                   Real zeta) const {
    assert(mmax <= mmax_);
    assert(T >= 0);
    const auto U = 0.25 * zeta * zeta * one_over_rho;
    assert(U >= 0);
    bool success = false;  // did interpolation or recursion work?
    if (T > Tmax || U < Umin) {
      success = eval_urr<false>(Gm, T, U, /* zeta_over_two_rho = */ 0, mmax);
    } else {
      success =
          interpolate_Gm<false>(Gm, T, U, /* zeta_over_two_rho = */ 0, mmax);
    }
    if (!success) {  // fallback is to use the Gaussian fit
      auto& gaussian_fit = get_gaussian_fit(zeta);
      if (!yukawa_gaussian_fit_eval_)
        yukawa_gaussian_fit_eval_ = std::unique_ptr<GaussianGmEval<Real, -1>>(
            new GaussianGmEval<Real, -1>(mmax_, precision_));
      yukawa_gaussian_fit_eval_->eval(Gm, 1.0 / one_over_rho, T, mmax,
                                      gaussian_fit);
    }
  }
  /// @param[in] Gm pointer to array of @c mmax+1 @c Real elements, on
  ///               return this contains the core integral for Slater-type
  ///               geminal, \f$ \exp(-zeta r_{12}) \f$ ,
  ///               namely \f$ (G_{m-1}(T,U) - G_m(T,U)) \zeta / (2 \rho), m \in
  ///               [0,m_\mathrm{max}] \f$ where \f$ G_{m}(T,U) = \int_0^1
  ///               t^{2m} \exp(U(1-t^{-2}) - Tt^2) dt \f$, where \f$ T = \rho
  ///               |\vec{P} - \vec{Q}|^2 \f$ and \f$ U = \zeta^2 / (4 \rho) \f$
  /// @param[in] one_over_rho \f$ 1/\rho \f$
  /// @param[in] T \f$ T \f$
  /// @param[in] mmax \f$ m_\mathrm{max} \f$
  /// @param[in] zeta \f$ \zeta \f$
  /// @throw std::invalid_argument if \f$ \zeta^2/(4 \rho) \f$ exceeds \c Umax
  void eval_slater(Real* Gm, Real one_over_rho, Real T, size_t mmax,
                   Real zeta) const {
    assert(mmax <= mmax_);
    assert(T >= 0);
    const auto U = 0.25 * zeta * zeta * one_over_rho;
    assert(U > 0);  // integral factorizes into 2 overlaps for U = 0
    const auto zeta_over_two_rho = 0.5 * zeta * one_over_rho;
    bool success = false;  // did interpolation or recursion work?
    if (T > Tmax) {
      success = eval_urr<true>(Gm, T, U, zeta_over_two_rho, mmax);
    } else {
      success = interpolate_Gm<true>(Gm, T, U, zeta_over_two_rho, mmax);
    }
    if (!success) {  // fallback is to use the Gaussian fit
      auto& gaussian_fit = get_gaussian_fit(zeta);
      if (!slater_gaussian_fit_eval_)
        slater_gaussian_fit_eval_ = std::unique_ptr<GaussianGmEval<Real, 0>>(
            new GaussianGmEval<Real, 0>(mmax_, precision_));
      slater_gaussian_fit_eval_->eval(Gm, 1.0 / one_over_rho, T, mmax,
                                      gaussian_fit);
    }
  }

  /// Fits a slater with an undetermined number of gaussians, borrowed from
  /// https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/madness/mra/gfit.h#L140
  /// @param[in] zeta exponent of Slater
  /// @param[in] lo smallest lengthscale to be represented
  /// @param[in] hi largest lengthscale to be represented
  /// @param[in] eps precision threshold
  /// @return vector of Gaussian geminals representing a Slater function with
  // the given zeta as a vector of std::pair(exponent, coefficient)
  static std::vector<std::pair<double, double>> fit_slater(double zeta,
                                                           double lo = 1e-5,
                                                           double hi = 5,
                                                           double eps = 1e-12) {
    double ulim;
    if (eps >= 1e-2)
      ulim = 5;
    else if (eps >= 1e-4)
      ulim = 10;
    else if (eps >= 1e-6)
      ulim = 14;
    else if (eps >= 1e-8)
      ulim = 18;
    else if (eps >= 1e-10)
      ulim = 22;
    else if (eps >= 1e-12)
      ulim = 26;
    else
      ulim = 30;

    // integration limits for quadrature
    double slo = 0.5 * log(eps) - 1.0;
    double shi = log(ulim / (lo * lo)) * 0.5;
    // resolution for quadrature
    double res = 1.0 / (0.2 - 0.5 * log10(eps));

    // Truncate the number of binary digits in h's mantissa
    // so that rounding does not occur when performing
    // manipulations to determine the quadrature points and
    // to limit the number of distinct values in case of
    // multiple precisions being used at the same time.
    res = floor(64. * res) / 64.;
    // round slo/hi down/up to an integral multiple of quadrature points
    slo = floor(slo / res) * res;
    shi = ceil(shi / res) * res;

    long n_pts = long((shi - slo) / res + 0.5);

    std::vector<double> coeffs(n_pts), expnts(n_pts);

    for (size_t i = 0; i != n_pts; ++i) {  // fill coeff and exponents vectors
      const double s = slo + res * (n_pts - i);
      coeffs[i] = res * exp(-zeta * zeta * exp(2. * s) + s);
      const long double twoOvrSqrtPi = 1.12837916709551257389615890312154517;
      coeffs[i] *= zeta * twoOvrSqrtPi;
      expnts[i] = 0.25 * exp(-2. * s);
    }

    assert(coeffs.size() == expnts.size());
    std::vector<std::pair<double, double>> result;
    for (size_t i = 0; i != coeffs.size(); ++i) {
      if (std::abs(coeffs[i]) > std::numeric_limits<double>::epsilon())
        result.push_back(std::pair<double, double>(expnts[i], coeffs[i]));
    }
    return result;
  }

 private:
  /// Computes G_0 and (optionally) G_{-1}
  /// @tparam compute_Gm1, if true, compute G_0 and G_{-1}, otherwise G_0 only
  /// @param[in] T the value of \f$ T \f$
  /// @param[in] U the value of \f$ U \f$
  /// @return if {T,U} can be handled without overflow return {G_0,G_{-1}} if
  //          @c compute_Gm1==true , else {G_0,0};
  //          if {T,U} cannot be handled without overflow
  //          return nullopt
  template <bool compute_Gm1>
  static inline optional<std::tuple<Real, Real>> eval_G0_and_maybe_Gm1(Real T,
                                                                       Real U) {
    assert(T >= 0);
    assert(U >= 0);
    const Real sqrtPi(1.7724538509055160272981674833411451827975494561224);

    Real G_m1 = 0;
    Real G_0 = 0;
    if (U ==
        0) {  // \sqrt{U} G_{-1} is finite, need to handle that case separately
      assert(!compute_Gm1);
      if (T < std::numeric_limits<Real>::epsilon()) {
        G_0 = 1;
      } else {
        const Real sqrt_T = sqrt(T);
        const Real sqrtPi_over_2 = sqrtPi * 0.5;
        G_0 = sqrtPi_over_2 * erf(sqrt_T) / sqrt_T;
      }
    } else if (T > std::numeric_limits<Real>::epsilon()) {  // U > 0
      const Real sqrt_U = sqrt(U);
      const Real sqrt_T = sqrt(T);
      const Real oo_sqrt_T = 1 / sqrt_T;
      const Real oo_sqrt_U = compute_Gm1 ? (1 / sqrt_U) : 0;
      const Real kappa = sqrt_U - sqrt_T;
      const Real lambda = sqrt_U + sqrt_T;
      const Real sqrtPi_over_4 = sqrtPi * 0.25;
      const Real pfac = sqrtPi_over_4;
      const Real erfc_k = exp(kappa * kappa - T) * erfc(kappa);
      // TODO when lambda * lambda - T exceeds 709 exp overflows, so need to
      // handle this more carefully
      if (lambda * lambda - T > 709) {
        libint2::verbose_stream()
            << "TennoGmEval::eval_G0_and_maybe_Gm1(): for large T and/or U "
            << "current implementation overflows, will represent "
            << "Slater function as a linear combination of Gaussians\n";
        return libint2::nullopt;
      }
      const Real erfc_l = exp(lambda * lambda - T) * erfc(lambda);

      G_m1 = compute_Gm1 ? pfac * (erfc_k + erfc_l) * oo_sqrt_U : 0.0;
      G_0 = pfac * (erfc_k - erfc_l) * oo_sqrt_T;
    } else {           // T = 0, U > 0
      if (U <= 709) {  // exp(<=709) is finite
        const Real exp_U = exp(U);
        const Real sqrt_U = sqrt(U);
        if (compute_Gm1) {
          const Real sqrtPi_over_2 = sqrtPi * 0.5;
          const Real oo_sqrt_U = compute_Gm1 ? (1 / sqrt_U) : 0;

          G_m1 = exp_U * sqrtPi_over_2 * erfc(sqrt_U) * oo_sqrt_U;
        }
        G_0 = 1 - exp_U * sqrtPi * erfc(sqrt_U) * sqrt_U;
      } else {  // exp(>=710) is "infinity", handle larger values as a series
        const Real oo_U = Real{1} / U;
        const Real x = oo_U;
        const Real x2 = x * x;
        const Real x3 = x2 * x;
        const Real x4 = x2 * x2;
        const Real x5 = x3 * x2;
        const Real x6 = x3 * x3;
        if (compute_Gm1) {
          // 6th-order is sufficient
          // in Wolfram: Series[Exp[U]*(Sqrt[Pi]/2)*Erfc[Sqrt[U]]/Sqrt[U], {U,
          // Infinity, 6}]
          const Real cm1[] = {1. / 2,    -1. / 4,   3. / 8,
                              -15. / 16, 105. / 32, -945. / 64};
          G_m1 = cm1[0] * x + cm1[1] * x2 + cm1[2] * x3 + cm1[3] * x4 +
                 cm1[4] * x5 + cm1[5] * x6;
        }
        // 6th-order is sufficient
        // in Wolfram: Series[1 - Exp[U]*Sqrt[Pi]*Erfc[Sqrt[U]]*Sqrt[U], {U,
        // Infinity, 6}]
        const Real c0[] = {1. / 2,     -3. / 4,   15. / 8,
                           -105. / 16, 945. / 32, -10395. / 64};
        G_0 = c0[0] * x + c0[1] * x2 + c0[2] * x3 + c0[3] * x4 + c0[4] * x5 +
              c0[5] * x6;
      }
    }

    return std::make_tuple(G_0, G_m1);
  }

  /// Computes core integrals by upward recursion from \f$ G_{-1} (T,U) \f$ and
  /// \f$ G_0 (T,U) \f$ this is unstable for small T, so use when interpolation
  /// does not apply.
  /// @tparam[in] Exp if true, will compute core integrals of the exponential
  /// operator, otherwise will compute Yukawa
  /// @param[out] Gm_vec pointer to an array with at least @c mmax+1 elements,
  /// on output contains:
  ///    - if `Exp==false`, `Gm_vec[m]` contains \f$ G_m(T,U), m=0..mmax\f$ ;
  ///    - if `Exp==true`, `Gm_vec[m]` contains \f$ \frac{\zeta}{2 \rho} \left(
  ///    G_{m-1}(T,U) - G_m(T,U) \right), m=0..mmax \f$;
  /// @param[in] T the value of \f$ T \f$
  /// @param[in] U the value of \f$ U \f$
  /// @param[in] zeta_over_two_rho the value of \f$ \frac{\zeta}{2 \rho} \f$,
  /// only used if @c Exp==true
  /// @param[in] mmax the maximum value of \f$ m \f$
  template <bool Exp>
  static bool eval_urr(Real* Gm_vec, Real T, Real U, Real zeta_over_two_rho,
                       size_t mmax) {
    assert(T >= 0);
    assert(U > 0);

    if (T == 0) return false;

    const Real sqrt_U = sqrt(U);
    const Real sqrt_T = sqrt(T);
    const Real oo_sqrt_T = 1 / sqrt_T;
    const Real oo_sqrt_U = 1 / sqrt_U;
    const Real kappa = sqrt_U - sqrt_T;
    const Real lambda = sqrt_U + sqrt_T;
    const Real sqrtPi_over_4(
        0.44311346272637900682454187083528629569938736403060);
    const Real pfac = sqrtPi_over_4;
    const Real erfc_k = exp(kappa * kappa - T) * erfc(kappa);
    const Real erfc_l = exp(lambda * lambda - T) * erfc(lambda);

    Real Gmm1;                                    // will contain G[m-1]
    Real Gm;                                      // will contain G[m]
    Real Gmp1;                                    // will contain G[m+1]
    Gmm1 = pfac * (erfc_k + erfc_l) * oo_sqrt_U;  // G_{-1}
    Gm = pfac * (erfc_k - erfc_l) * oo_sqrt_T;    // G_{0}
    if
#if __cplusplus >= 201703L
        constexpr
#endif
        (!Exp) {
      Gm_vec[0] = Gm;
    } else {
      Gm_vec[0] = (Gmm1 - Gm) * zeta_over_two_rho;
    }

    if (mmax > 0) {
      // first application of URR
      const Real oo_two_T = 0.5 / T;
      const Real two_U = 2.0 * U;
      const Real exp_mT = exp(-T);

      for (unsigned int m = 0, two_m_minus_1 = 1; m < mmax;
           ++m, two_m_minus_1 += 2) {
        Gmp1 = oo_two_T * (two_m_minus_1 * Gm + two_U * Gmm1 - exp_mT);
        if
#if __cplusplus >= 201703L
            constexpr
#endif
            (!Exp) {
          Gm_vec[m + 1] = Gmp1;
        } else {
          Gm_vec[m + 1] = (Gm - Gmp1) * zeta_over_two_rho;
        }
        Gmm1 = Gm;
        Gm = Gmp1;
      }
    }

    return true;
  }

  /// compute core integrals by Chebyshev interpolation
  /// @tparam[in] Exp if true, will compute core integrals of the exponential
  /// operator, otherwise will compute Yukawa
  /// @param[out] Gm_vec if @c Exp==false Gm_vec[m]=\f$ G_m(T,U), m=0..mmax \f$,
  /// otherwise Gm_vec[m]=\f$ \frac{\zeta}{2 \rho} \left( G_{m-1}(T,U) -
  /// G_m(T,U) \right), m=0..mmax \f$ ; must be at least @c mmax+1 elements long
  /// @param[in] T the value of \f$ T \f$
  /// @param[in] U the value of \f$ U \f$
  /// @param[in] zeta_over_two_rho the value of \f$ \frac{\zeta}{2 \rho} \f$,
  /// only used if @c Exp==true
  /// @param[in] mmax the maximum value of \f$ m \f$
  /// @throw std::invalid_argument if \p U exceeds \c Umax
  /// @return true if {T,U} can be handled without overflow; false return
  ///         indicates that Gaussian fit of the kernel should be used to
  ///         compute Gm
  template <bool Exp>
  bool interpolate_Gm(Real* Gm_vec, Real T, Real U, Real zeta_over_two_rho,
                      long mmax) const {
    assert(T >= 0);
    if (U > Umax) {
      libint2::verbose_stream()
          << "TennoGmEval::eval_{yukawa,slater}() : arguments out of range, "
          << "zeta*zeta*one_over_rho/4=" << std::to_string(U)
          << " cannot exceed " << std::to_string(Umax)
          << " will represent Slater function as a linear combination of"
          << " Gaussian functions " << std::endl;
      return false;
    }
    if (U < Umin || U > Umax) {
      return false;
    }

    // maps x in [0,2] to [-1/2,1/2] in a linear fashion
    auto linear_map_02 = [](Real x) {
      assert(x >= 0 && x <= 2);
      return (x - 1) * 0.5;
    };
    // maps x in [a, 2*a] to [-1/2,1/2] in a logarithmic fashion
    auto log2_map = [](Real x, Real one_over_a) {
      return std::log2(x * one_over_a) - 0.5;
    };
    // maps x in [a, 10*a] to [-1/2,1/2] in a logarithmic fashion
    auto log10_map = [](Real x, Real one_over_a) {
      return std::log10(x * one_over_a) - 0.5;
    };

    // which interval does this T fall into?
    const int Tint =
        T < 2 ? 0
              : int(std::floor(std::log2(T) *
                               (1. - std::numeric_limits<double>::epsilon())));
    assert(Tint >= 0 && Tint < cheb_table_tmaxlog2 - cheb_table_tminlog2 + 1);
    // precomputed 1 / 2^K
    constexpr Real one_over_2K[] = {1,          0.5,       0.25,     0.125,
                                    0.0625,     0.03125,   0.015625, 0.0078125,
                                    0.00390625, .001953125};
    // which interval does this U fall into?
    const auto log10_U_over_Umin =
        std::log10(U * one_over_Umin) *
        (1. - std::numeric_limits<double>::epsilon());
    const int Uint = int(std::floor(log10_U_over_Umin));
    assert(Uint >= 0 && Uint < cheb_table_umaxlog10 - cheb_table_uminlog10);
    // precomputed 1 / 10^(cheb_table_uminlog10 + K)
    constexpr Real one_over_10K[] = {
        1. / detail::pow10(cheb_table_uminlog10),
        1. / detail::pow10(cheb_table_uminlog10 + 1),
        1. / detail::pow10(cheb_table_uminlog10 + 2),
        1. / detail::pow10(cheb_table_uminlog10 + 3),
        1. / detail::pow10(cheb_table_uminlog10 + 4),
        1. / detail::pow10(cheb_table_uminlog10 + 5),
        1. / detail::pow10(cheb_table_uminlog10 + 6),
        1. / detail::pow10(cheb_table_uminlog10 + 7),
        1. / detail::pow10(cheb_table_uminlog10 + 8),
        1. / detail::pow10(cheb_table_uminlog10 + 9)};
    const Real t =
        Tint == 0
            ? linear_map_02(T)
            : log2_map(T, one_over_2K[Tint]);  // this ranges from -0.5 to 0.5
    const Real u =
        log10_map(U, one_over_10K[Uint]);  // this ranges from -0.5 to 0.5

    const int interval =
        Tint * (cheb_table_umaxlog10 - cheb_table_uminlog10) + Uint;

#if defined(__AVX__)
    assert(ORDERp1 == 16);
    const auto t2 = t * t;
    const auto t3 = t2 * t;
    const auto t4 = t2 * t2;
    const auto t5 = t2 * t3;
    const auto t6 = t3 * t3;
    const auto t7 = t3 * t4;
    const auto t8 = t4 * t4;
    const auto t9 = t4 * t5;
    const auto t10 = t5 * t5;
    const auto t11 = t5 * t6;
    const auto t12 = t6 * t6;
    const auto t13 = t6 * t7;
    const auto t14 = t7 * t7;
    const auto t15 = t7 * t8;
    const auto u2 = u * u;
    const auto u3 = u2 * u;
    const auto u4 = u2 * u2;
    const auto u5 = u2 * u3;
    const auto u6 = u3 * u3;
    const auto u7 = u3 * u4;
    const auto u8 = u4 * u4;
    const auto u9 = u4 * u5;
    const auto u10 = u5 * u5;
    const auto u11 = u5 * u6;
    const auto u12 = u6 * u6;
    const auto u13 = u6 * u7;
    const auto u14 = u7 * u7;
    const auto u15 = u7 * u8;
    libint2::simd::VectorAVXDouble u0vec(1., u, u2, u3);
    libint2::simd::VectorAVXDouble u1vec(u4, u5, u6, u7);
    libint2::simd::VectorAVXDouble u2vec(u8, u9, u10, u11);
    libint2::simd::VectorAVXDouble u3vec(u12, u13, u14, u15);
#endif  // AVX

    Real Gmm1 = 0.0;  // will track the previous value, only used for Exp==false

    constexpr bool compute_Gmm10 = true;
    long mmin_interp;
    if (compute_Gmm10) {
      // precision of interpolation for m=-1,0 can be insufficient, just
      // evaluate explicitly
      Real G0;
      const auto eval_G0_Gm1_result = eval_G0_and_maybe_Gm1<Exp>(T, U);
      if (!eval_G0_Gm1_result.has_value()) {
        return false;
      }
      std::tie(Exp ? G0 : Gm_vec[0], Gmm1) = eval_G0_Gm1_result.value();
      if
#if __cplusplus >= 201703L
          constexpr
#endif
          (Exp) {
        Gm_vec[0] = (Gmm1 - G0) * zeta_over_two_rho;
        Gmm1 = G0;
      }
      mmin_interp = 1;
    } else
      mmin_interp = Exp ? -1 : 0;

    // now compute the rest
    for (long m = mmin_interp; m <= mmax; ++m) {
      const Real* c_tuint =
          c_ + (ORDERp1) * (ORDERp1) *
                   (interval * (mmax_ - mmin_ + 1) +
                    (m - mmin_));  // ptr to the interpolation data for m=mmin
#if defined(__AVX__)
      libint2::simd::VectorAVXDouble c00v, c01v, c02v, c03v, c10v, c11v, c12v,
          c13v, c20v, c21v, c22v, c23v, c30v, c31v, c32v, c33v, c40v, c41v,
          c42v, c43v, c50v, c51v, c52v, c53v, c60v, c61v, c62v, c63v, c70v,
          c71v, c72v, c73v;
      libint2::simd::VectorAVXDouble c80v, c81v, c82v, c83v, c90v, c91v, c92v,
          c93v, ca0v, ca1v, ca2v, ca3v, cb0v, cb1v, cb2v, cb3v, cc0v, cc1v,
          cc2v, cc3v, cd0v, cd1v, cd2v, cd3v, ce0v, ce1v, ce2v, ce3v, cf0v,
          cf1v, cf2v, cf3v;
      c00v.load_aligned(c_tuint);
      c01v.load_aligned((c_tuint + 4));
      c02v.load_aligned(c_tuint + 8);
      c03v.load_aligned((c_tuint + 12));
      libint2::simd::VectorAVXDouble t0vec(1, 1, 1, 1);
      libint2::simd::VectorAVXDouble t0 =
          t0vec * (c00v * u0vec + c01v * u1vec + c02v * u2vec + c03v * u3vec);
      c10v.load_aligned(c_tuint + ORDERp1);
      c11v.load_aligned((c_tuint + 4) + ORDERp1);
      c12v.load_aligned((c_tuint + 8) + ORDERp1);
      c13v.load_aligned((c_tuint + 12) + ORDERp1);
      libint2::simd::VectorAVXDouble t1vec(t, t, t, t);
      libint2::simd::VectorAVXDouble t1 =
          t1vec * (c10v * u0vec + c11v * u1vec + c12v * u2vec + c13v * u3vec);
      c20v.load_aligned(c_tuint + 2 * ORDERp1);
      c21v.load_aligned((c_tuint + 4) + 2 * ORDERp1);
      c22v.load_aligned((c_tuint + 8) + 2 * ORDERp1);
      c23v.load_aligned((c_tuint + 12) + 2 * ORDERp1);
      libint2::simd::VectorAVXDouble t2vec(t2, t2, t2, t2);
      libint2::simd::VectorAVXDouble t2 =
          t2vec * (c20v * u0vec + c21v * u1vec + c22v * u2vec + c23v * u3vec);
      c30v.load_aligned(c_tuint + 3 * ORDERp1);
      c31v.load_aligned((c_tuint + 4) + 3 * ORDERp1);
      c32v.load_aligned((c_tuint + 8) + 3 * ORDERp1);
      c33v.load_aligned((c_tuint + 12) + 3 * ORDERp1);
      libint2::simd::VectorAVXDouble t3vec(t3, t3, t3, t3);
      libint2::simd::VectorAVXDouble t3 =
          t3vec * (c30v * u0vec + c31v * u1vec + c32v * u2vec + c33v * u3vec);
      libint2::simd::VectorAVXDouble t0123 = horizontal_add(t0, t1, t2, t3);

      c40v.load_aligned(c_tuint + 4 * ORDERp1);
      c41v.load_aligned((c_tuint + 4) + 4 * ORDERp1);
      c42v.load_aligned((c_tuint + 8) + 4 * ORDERp1);
      c43v.load_aligned((c_tuint + 12) + 4 * ORDERp1);
      libint2::simd::VectorAVXDouble t4vec(t4, t4, t4, t4);
      libint2::simd::VectorAVXDouble t4 =
          t4vec * (c40v * u0vec + c41v * u1vec + c42v * u2vec + c43v * u3vec);
      c50v.load_aligned(c_tuint + 5 * ORDERp1);
      c51v.load_aligned((c_tuint + 4) + 5 * ORDERp1);
      c52v.load_aligned((c_tuint + 8) + 5 * ORDERp1);
      c53v.load_aligned((c_tuint + 12) + 5 * ORDERp1);
      libint2::simd::VectorAVXDouble t5vec(t5, t5, t5, t5);
      libint2::simd::VectorAVXDouble t5 =
          t5vec * (c50v * u0vec + c51v * u1vec + c52v * u2vec + c53v * u3vec);
      c60v.load_aligned(c_tuint + 6 * ORDERp1);
      c61v.load_aligned((c_tuint + 4) + 6 * ORDERp1);
      c62v.load_aligned((c_tuint + 8) + 6 * ORDERp1);
      c63v.load_aligned((c_tuint + 12) + 6 * ORDERp1);
      libint2::simd::VectorAVXDouble t6vec(t6, t6, t6, t6);
      libint2::simd::VectorAVXDouble t6 =
          t6vec * (c60v * u0vec + c61v * u1vec + c62v * u2vec + c63v * u3vec);
      c70v.load_aligned(c_tuint + 7 * ORDERp1);
      c71v.load_aligned((c_tuint + 4) + 7 * ORDERp1);
      c72v.load_aligned((c_tuint + 8) + 7 * ORDERp1);
      c73v.load_aligned((c_tuint + 12) + 7 * ORDERp1);
      libint2::simd::VectorAVXDouble t7vec(t7, t7, t7, t7);
      libint2::simd::VectorAVXDouble t7 =
          t7vec * (c70v * u0vec + c71v * u1vec + c72v * u2vec + c73v * u3vec);
      libint2::simd::VectorAVXDouble t4567 = horizontal_add(t4, t5, t6, t7);

      c80v.load_aligned(c_tuint + 8 * ORDERp1);
      c81v.load_aligned((c_tuint + 4) + 8 * ORDERp1);
      c82v.load_aligned((c_tuint + 8) + 8 * ORDERp1);
      c83v.load_aligned((c_tuint + 12) + 8 * ORDERp1);
      libint2::simd::VectorAVXDouble t8vec(t8, t8, t8, t8);
      libint2::simd::VectorAVXDouble t8 =
          t8vec * (c80v * u0vec + c81v * u1vec + c82v * u2vec + c83v * u3vec);
      c90v.load_aligned(c_tuint + 9 * ORDERp1);
      c91v.load_aligned((c_tuint + 4) + 9 * ORDERp1);
      c92v.load_aligned((c_tuint + 8) + 9 * ORDERp1);
      c93v.load_aligned((c_tuint + 12) + 9 * ORDERp1);
      libint2::simd::VectorAVXDouble t9vec(t9, t9, t9, t9);
      libint2::simd::VectorAVXDouble t9 =
          t9vec * (c90v * u0vec + c91v * u1vec + c92v * u2vec + c93v * u3vec);
      ca0v.load_aligned(c_tuint + 10 * ORDERp1);
      ca1v.load_aligned((c_tuint + 4) + 10 * ORDERp1);
      ca2v.load_aligned((c_tuint + 8) + 10 * ORDERp1);
      ca3v.load_aligned((c_tuint + 12) + 10 * ORDERp1);
      libint2::simd::VectorAVXDouble tavec(t10, t10, t10, t10);
      libint2::simd::VectorAVXDouble ta =
          tavec * (ca0v * u0vec + ca1v * u1vec + ca2v * u2vec + ca3v * u3vec);
      cb0v.load_aligned(c_tuint + 11 * ORDERp1);
      cb1v.load_aligned((c_tuint + 4) + 11 * ORDERp1);
      cb2v.load_aligned((c_tuint + 8) + 11 * ORDERp1);
      cb3v.load_aligned((c_tuint + 12) + 11 * ORDERp1);
      libint2::simd::VectorAVXDouble tbvec(t11, t11, t11, t11);
      libint2::simd::VectorAVXDouble tb =
          tbvec * (cb0v * u0vec + cb1v * u1vec + cb2v * u2vec + cb3v * u3vec);
      libint2::simd::VectorAVXDouble t89ab = horizontal_add(t8, t9, ta, tb);

      cc0v.load_aligned(c_tuint + 12 * ORDERp1);
      cc1v.load_aligned((c_tuint + 4) + 12 * ORDERp1);
      cc2v.load_aligned((c_tuint + 8) + 12 * ORDERp1);
      cc3v.load_aligned((c_tuint + 12) + 12 * ORDERp1);
      libint2::simd::VectorAVXDouble tcvec(t12, t12, t12, t12);
      libint2::simd::VectorAVXDouble tc =
          tcvec * (cc0v * u0vec + cc1v * u1vec + cc2v * u2vec + cc3v * u3vec);
      cd0v.load_aligned(c_tuint + 13 * ORDERp1);
      cd1v.load_aligned((c_tuint + 4) + 13 * ORDERp1);
      cd2v.load_aligned((c_tuint + 8) + 13 * ORDERp1);
      cd3v.load_aligned((c_tuint + 12) + 13 * ORDERp1);
      libint2::simd::VectorAVXDouble tdvec(t13, t13, t13, t13);
      libint2::simd::VectorAVXDouble td =
          tdvec * (cd0v * u0vec + cd1v * u1vec + cd2v * u2vec + cd3v * u3vec);
      ce0v.load_aligned(c_tuint + 14 * ORDERp1);
      ce1v.load_aligned((c_tuint + 4) + 14 * ORDERp1);
      ce2v.load_aligned((c_tuint + 8) + 14 * ORDERp1);
      ce3v.load_aligned((c_tuint + 12) + 14 * ORDERp1);
      libint2::simd::VectorAVXDouble tevec(t14, t14, t14, t14);
      libint2::simd::VectorAVXDouble te =
          tevec * (ce0v * u0vec + ce1v * u1vec + ce2v * u2vec + ce3v * u3vec);
      cf0v.load_aligned(c_tuint + 15 * ORDERp1);
      cf1v.load_aligned((c_tuint + 4) + 15 * ORDERp1);
      cf2v.load_aligned((c_tuint + 8) + 15 * ORDERp1);
      cf3v.load_aligned((c_tuint + 4) + 15 * ORDERp1);
      libint2::simd::VectorAVXDouble tfvec(t15, t15, t15, t15);
      libint2::simd::VectorAVXDouble tf =
          tfvec * (cf0v * u0vec + cf1v * u1vec + cf2v * u2vec + cf3v * u3vec);
      libint2::simd::VectorAVXDouble tcdef = horizontal_add(tc, td, te, tf);

      auto tall = horizontal_add(t0123, t4567, t89ab, tcdef);
      const auto Gm = horizontal_add(tall);
#else   // AVX
      // no AVX, use explicit loops (probably slow)
      double uvec[16];
      double tvec[16];
      uvec[0] = 1;
      tvec[0] = 1;
      for (int i = 1; i != 16; ++i) {
        uvec[i] = uvec[i - 1] * u;
        tvec[i] = tvec[i - 1] * t;
      }
      double Gm = 0.0;
      for (int i = 0, ij = 0; i != 16; ++i) {
        for (int j = 0; j != 16; ++j, ++ij) {
          Gm += c_tuint[ij] * tvec[i] * uvec[j];
        }
      }
#endif  // AVX

      if
#if __cplusplus >= 201703L
          constexpr
#endif
          (!Exp) {
        Gm_vec[m] = Gm;
      } else {
        if (m != -1) {
          Gm_vec[m] = (Gmm1 - Gm) * zeta_over_two_rho;
        }
        Gmm1 = Gm;
      }
    }

    return true;
  }

 private:
  unsigned int mmax_;
  Real precision_;
  ExpensiveNumbers<Real> numbers_;
  Real* c_ = nullptr;

  void init_table() {
    // get memory
    void* result;
    int status = posix_memalign(&result, std::max(sizeof(Real), (size_t)32),
                                (mmax_ - mmin_ + 1) * cheb_table_nintervals *
                                    ORDERp1 * ORDERp1 * sizeof(Real));
    if (status != 0) {
      if (status == EINVAL)
        throw std::logic_error(
            "TennoGmEval::init() : posix_memalign failed, alignment must be "
            "a "
            "power of 2 at least as large as sizeof(void *)");
      if (status == ENOMEM) throw std::bad_alloc();
      abort();  // should be unreachable
    }
    c_ = static_cast<Real*>(result);

    // copy contents of static table into c
    // need all intervals
    for (std::size_t iv = 0; iv < cheb_table_nintervals; ++iv) {
      // but only values of m up to mmax
      std::copy(cheb_table[iv],
                cheb_table[iv] + ((mmax_ - mmin_) + 1) * ORDERp1 * ORDERp1,
                c_ + (iv * ((mmax_ - mmin_) + 1)) * ORDERp1 * ORDERp1);
    }
  }

};  // TennoGmEval

#if LIBINT2_CONSTEXPR_STATICS
template <typename Real>
constexpr double
    TennoGmEval<Real>::cheb_table[TennoGmEval<Real>::cheb_table_nintervals]
                                 [(TennoGmEval<Real>::cheb_table_mmax + 2) *
                                  (TennoGmEval<Real>::interpolation_order + 1) *
                                  (TennoGmEval<Real>::interpolation_order + 1)];
#else
// clang needs an explicit specialization declaration to avoid warning
#ifdef __clang__
template <typename Real>
double
    TennoGmEval<Real>::cheb_table[TennoGmEval<Real>::cheb_table_nintervals]
                                 [(TennoGmEval<Real>::cheb_table_mmax + 2) *
                                  (TennoGmEval<Real>::interpolation_order + 1) *
                                  (TennoGmEval<Real>::interpolation_order + 1)];
#endif
#endif

template <typename Real, int k>
struct GaussianGmEval;

namespace detail {

/// some evaluators need thread-local scratch, but most don't
template <typename CoreEval>
struct CoreEvalScratch {
  CoreEvalScratch(const CoreEvalScratch&) = default;
  CoreEvalScratch(CoreEvalScratch&&) = default;
  explicit CoreEvalScratch(int) {}
};
/// GaussianGmEval<Real,-1> needs extra scratch data
template <typename Real>
struct CoreEvalScratch<GaussianGmEval<Real, -1>> {
  std::vector<Real> Fm_;
  std::vector<Real> g_i;  // (gamma / (rho  + gamma))^i
  std::vector<Real> r_i;  // (rho / (rho  + gamma))^i
  CoreEvalScratch(const CoreEvalScratch&) = default;
  CoreEvalScratch(CoreEvalScratch&&) = default;
  explicit CoreEvalScratch(int mmax) { init(mmax); }

 private:
  void init(int mmax) {
    Fm_.resize(mmax + 1);
    g_i.resize(mmax + 1);
    r_i.resize(mmax + 1);
    g_i[0] = 1.0;
    r_i[0] = 1.0;
  }
};
}  // namespace detail

//////////////////////////////////////////////////////////
/// core integrals r12^k \sum_i \exp(- a_i r_12^2)
//////////////////////////////////////////////////////////

/**
 * Evaluates core integral \$ G_m(\rho, T) = \left( - \frac{\partial}{\partial
 * T} \right)^n G_0(\rho,T) \f$, \f$ G_0(\rho,T) = \int \exp(-\rho |\vec{r} -
 * \vec{P} + \vec{Q}|^2) g(r) \, {\rm d}\vec{r} \f$ over a general contracted
 * Gaussian geminal \f$ g(r_{12}) = r_{12}^k \sum_i c_i \exp(- a_i r_{12}^2),
 * \quad k = -1, 0, 2 \f$ . The integrals are needed in R12/F12 methods with
 * STG-nG correlation factors. Specifically, for a correlation factor \f$
 * f(r_{12}) = \sum_i c_i \exp(- a_i r_{12}^2) \f$ integrals with the
 * following kernels are needed: <ul> <li> \f$ f(r_{12}) \f$  (k=0) </li> <li>
 * \f$ f(r_{12}) / r_{12} \f$  (k=-1) </li> <li> \f$ f(r_{12})^2 \f$ (k=0, @sa
 * GaussianGmEval::eval ) </li> <li> \f$ [f(r_{12}), [\hat{T}_1, f(r_{12})]]
 * \f$ (k=2, @sa GaussianGmEval::eval ) </li>
 * </ul>
 *
 * N.B. ``Asymmetric'' kernels, \f$ f(r_{12}) g(r_{12}) \f$ and
 *   \f$ [f(r_{12}), [\hat{T}_1, g(r_{12})]] \f$, where f and g are two
 * different geminals, can also be handled straightforwardly.
 *
 * \note for more details see DOI: 10.1039/b605188j
 */
template <typename Real, int k>
struct GaussianGmEval
    : private detail::CoreEvalScratch<
          GaussianGmEval<Real, k>>  // N.B. empty-base optimization
{
#ifndef LIBINT_USER_DEFINED_REAL
  using FmEvalType = libint2::FmEval_Chebyshev7<double>;
#else
  using FmEvalType = libint2::FmEval_Reference<scalar_type>;
#endif

  /**
   * @param[in] mmax the evaluator will be used to compute Gm(T) for 0 <= m <=
   * mmax
   */
  GaussianGmEval(int mmax, Real precision)
      : detail::CoreEvalScratch<GaussianGmEval<Real, k>>(mmax),
        mmax_(mmax),
        precision_(precision),
        fm_eval_(),
        numbers_(-1, -1, mmax) {
    assert(k == -1 || k == 0 || k == 2);
    // for k=-1 need to evaluate the Boys function
    if (k == -1) {
      fm_eval_ = FmEvalType::instance(mmax_, precision_);
    }
  }

  ~GaussianGmEval() {}

  /// Singleton interface allows to manage the lone instance;
  /// adjusts max m and precision values as needed in thread-safe fashion
  static std::shared_ptr<GaussianGmEval> instance(
      unsigned int mmax,
      Real precision = std::numeric_limits<Real>::epsilon()) {
    assert(mmax >= 0);
    assert(precision >= 0);
    // thread-safe per C++11 standard [6.7.4]
    static auto instance_ = std::make_shared<GaussianGmEval>(mmax, precision);

    while (instance_->max_m() < mmax || instance_->precision() > precision) {
      static std::mutex mtx;
      std::lock_guard<std::mutex> lck(mtx);
      if (instance_->max_m() < mmax || instance_->precision() > precision) {
        auto new_instance = std::make_shared<GaussianGmEval>(mmax, precision);
        instance_ = new_instance;
      }
    }

    return instance_;
  }

  /// @return the maximum value of m for which the \f$ G_m(\rho, T) \f$ can be
  /// computed with this object
  int max_m() const { return mmax_; }
  /// @return the precision with which this object can compute the Boys
  /// function
  Real precision() const { return precision_; }

  /** computes \f$ G_m(\rho, T) \f$ using downward recursion.
   *
   * @warning NOT reentrant if \c k == -1 and C++11 is not available
   *
   * @param[out] Gm array to be filled in with the \f$ Gm(\rho, T) \f$ values,
   * must be at least mmax+1 elements long
   * @param[in] rho
   * @param[in] T
   * @param[in] mmax mmax the maximum value of m for which Boys function will
   * be computed; it must be <= the value returned by max_m() (this is not
   * checked)
   * @param[in] geminal the Gaussian geminal for which the core integral \f$
   * Gm(\rho, T) \f$ is computed; a contracted Gaussian geminal is represented
   * as a vector of pairs, each of which specifies the exponent (first) and
   * contraction coefficient (second) of the primitive Gaussian geminal
   * @param[in] scr if \c k ==-1 and need this to be reentrant, must provide
   * ptr to the per-thread \c
   * libint2::detail::CoreEvalScratch<GaussianGmEval<Real,-1>> object; no need
   * to specify \c scr otherwise
   */
  template <typename AnyReal>
  void eval(Real* Gm, Real rho, Real T, size_t mmax,
            const std::vector<std::pair<AnyReal, AnyReal>>& geminal,
            void* scr = 0) {
    std::fill(Gm, Gm + mmax + 1, Real(0));

    const auto sqrt_rho = sqrt(rho);
    const auto oo_sqrt_rho = 1 / sqrt_rho;

    typedef
        typename std::vector<std::pair<AnyReal, AnyReal>>::const_iterator citer;
    const citer gend = geminal.end();
    for (citer i = geminal.begin(); i != gend; ++i) {
      const auto gamma = i->first;
      const auto gcoef = i->second;
      const auto rhog = rho + gamma;
      const auto oorhog = 1 / rhog;

      const auto gorg = gamma * oorhog;
      const auto rorg = rho * oorhog;
      const auto sqrt_rho_org = sqrt_rho * oorhog;
      const auto sqrt_rhog = sqrt(rhog);
      const auto sqrt_rorg = sqrt_rho_org * sqrt_rhog;

      /// (ss|g12|ss)
      constexpr Real const_SQRTPI_2(
          0.88622692545275801364908374167057259139877472806119); /* sqrt(pi)/2
                                                                  */
      const auto SS_K0G12_SS = gcoef * oo_sqrt_rho * const_SQRTPI_2 * rorg *
                               sqrt_rorg * exp(-gorg * T);

      if (k == -1) {
        void* _scr = (scr == 0) ? this : scr;
        auto& scratch =
            *(reinterpret_cast<
                detail::CoreEvalScratch<GaussianGmEval<Real, -1>>*>(_scr));

        const auto rorgT = rorg * T;
        fm_eval_->eval(&scratch.Fm_[0], rorgT, mmax);

#if 1
        constexpr Real const_2_SQRTPI(
            1.12837916709551257389615890312154517); /* 2/sqrt(pi)     */
        const auto pfac = const_2_SQRTPI * sqrt_rhog * SS_K0G12_SS;
        for (int i = 1; i <= mmax; i++) {
          scratch.r_i[i] = scratch.r_i[i - 1] * rorg;
          scratch.g_i[i] = scratch.g_i[i - 1] * gorg;
        }
        for (int m = 0; m <= mmax; m++) {
          Real ssss = 0.0;
          Real* bcm = numbers_.bc[m];
          for (int n = 0; n <= m; n++) {
            ssss +=
                bcm[n] * scratch.r_i[n] * scratch.g_i[m - n] * scratch.Fm_[n];
          }
          Gm[m] += ssss * pfac;
        }
#endif
      }

      if (k == 0) {
        auto ss_oper_ss_m = SS_K0G12_SS;
        Gm[0] += ss_oper_ss_m;
        for (int m = 1; m <= mmax; ++m) {
          ss_oper_ss_m *= gorg;
          Gm[m] += ss_oper_ss_m;
        }
      }

      if (k == 2) {
        /// (ss|g12*r12^2|ss)
        const auto rorgT = rorg * T;
        const auto SS_K2G12_SS_0 = (1.5 + rorgT) * (SS_K0G12_SS * oorhog);
        const auto SS_K2G12_SS_m1 = rorg * (SS_K0G12_SS * oorhog);

        auto SS_K2G12_SS_gorg_m = SS_K2G12_SS_0;
        auto SS_K2G12_SS_gorg_m1 = SS_K2G12_SS_m1;
        Gm[0] += SS_K2G12_SS_gorg_m;
        for (int m = 1; m <= mmax; ++m) {
          SS_K2G12_SS_gorg_m *= gorg;
          Gm[m] += SS_K2G12_SS_gorg_m - m * SS_K2G12_SS_gorg_m1;
          SS_K2G12_SS_gorg_m1 *= gorg;
        }
      }
    }
  }

 private:
  int mmax_;
  Real precision_;  //< absolute precision
  std::shared_ptr<const FmEvalType> fm_eval_;

  ExpensiveNumbers<Real> numbers_;
};

template <typename GmEvalFunction>
struct GenericGmEval : private GmEvalFunction {
  typedef typename GmEvalFunction::value_type Real;

  GenericGmEval(int mmax, Real precision)
      : GmEvalFunction(mmax, precision), mmax_(mmax), precision_(precision) {}

  static std::shared_ptr<const GenericGmEval> instance(
      int mmax, Real precision = std::numeric_limits<Real>::epsilon()) {
    return std::make_shared<const GenericGmEval>(mmax, precision);
  }

  template <typename Real, typename... ExtraArgs>
  void eval(Real* Gm, Real rho, Real T, int mmax, ExtraArgs... args) const {
    assert(mmax <= mmax_);
    (GmEvalFunction(*this))(Gm, rho, T, mmax, std::forward<ExtraArgs>(args)...);
  }

  /// @return the maximum value of m for which the \f$ G_m(\rho, T) \f$ can be
  /// computed with this object
  int max_m() const { return mmax_; }
  /// @return the precision with which this object can compute the Boys
  /// function
  Real precision() const { return precision_; }

 private:
  int mmax_;
  Real precision_;
};

// these Gm engines need extra scratch data
namespace os_core_ints {
template <typename Real, int K>
struct r12_xx_K_gm_eval;
template <typename Real>
struct erfx_coulomb_gm_eval;
}  // namespace os_core_ints

namespace detail {
/// r12_xx_K_gm_eval<1> needs extra scratch data
template <typename Real>
struct CoreEvalScratch<os_core_ints::r12_xx_K_gm_eval<Real, 1>> {
  std::vector<Real> Fm_;
  CoreEvalScratch(const CoreEvalScratch&) = default;
  CoreEvalScratch(CoreEvalScratch&&) = default;
  // need to store Fm(T) for m = 0 .. mmax+1
  explicit CoreEvalScratch(int mmax) { Fm_.resize(mmax + 2); }
};
/// erfx_coulomb_gm_eval needs extra scratch data
template <typename Real>
struct CoreEvalScratch<os_core_ints::erfx_coulomb_gm_eval<Real>> {
  std::vector<Real> Fm_;
  CoreEvalScratch(const CoreEvalScratch&) = default;
  CoreEvalScratch(CoreEvalScratch&&) = default;
  // need to store Fm(T) for m = 0 .. mmax
  explicit CoreEvalScratch(int mmax) { Fm_.resize(mmax + 1); }
};
}  // namespace detail

/// Obara-Saika core ints code
namespace os_core_ints {

/// core integral evaluator delta function kernels
template <typename Real>
struct delta_gm_eval {
  typedef Real value_type;

  delta_gm_eval(unsigned int, Real) {}
  void operator()(Real* Gm, Real rho, Real T, int mmax) const {
    constexpr static auto one_over_two_pi = 1.0 / (2.0 * M_PI);
    const auto G0 = exp(-T) * rho * one_over_two_pi;
    std::fill(Gm, Gm + mmax + 1, G0);
  }
};

/// core integral evaluator for \f$ r_{12}^K \f$ kernel
/// @tparam K currently supported \c K=1 (use Boys engine directly for \c
/// K=-1)
/// @note need extra scratch for Boys function values when \c K==1,
///       the Gm vector is not long enough for scratch

template <typename Real, int K>
struct r12_xx_K_gm_eval;

template <typename Real>
struct r12_xx_K_gm_eval<Real, 1>
    : private detail::CoreEvalScratch<r12_xx_K_gm_eval<Real, 1>> {
  typedef detail::CoreEvalScratch<r12_xx_K_gm_eval<Real, 1>> base_type;
  typedef Real value_type;

#ifndef LIBINT_USER_DEFINED_REAL
  using FmEvalType = libint2::FmEval_Chebyshev7<double>;
#else
  using FmEvalType = libint2::FmEval_Reference<scalar_type>;
#endif

  r12_xx_K_gm_eval(unsigned int mmax, Real precision) : base_type(mmax) {
    fm_eval_ = FmEvalType::instance(mmax + 1, precision);
  }
  void operator()(Real* Gm, Real rho, Real T, int mmax) {
    fm_eval_->eval(&base_type::Fm_[0], T, mmax + 1);
    auto T_plus_m_plus_one = T + 1.0;
    Gm[0] = T_plus_m_plus_one * base_type::Fm_[0] - T * base_type::Fm_[1];
    auto minus_m = -1.0;
    T_plus_m_plus_one += 1.0;
    for (auto m = 1; m <= mmax; ++m, minus_m -= 1.0, T_plus_m_plus_one += 1.0) {
      Gm[m] = minus_m * base_type::Fm_[m - 1] +
              T_plus_m_plus_one * base_type::Fm_[m] - T * base_type::Fm_[m + 1];
    }
  }

 private:
  std::shared_ptr<const FmEvalType> fm_eval_;  // need for odd K
};

/// core integral evaluator for \f$ \mathrm{erf}(\omega r) / r \f$ kernel
template <typename Real>
struct erf_coulomb_gm_eval {
  typedef Real value_type;

#ifndef LIBINT_USER_DEFINED_REAL
  using FmEvalType = libint2::FmEval_Chebyshev7<double>;
#else
  using FmEvalType = libint2::FmEval_Reference<scalar_type>;
#endif

  erf_coulomb_gm_eval(unsigned int mmax, Real precision) {
    fm_eval_ = FmEvalType::instance(mmax, precision);
  }
  void operator()(Real* Gm, Real rho, Real T, int mmax, Real omega) const {
    if (omega > 0) {
      auto omega2 = omega * omega;
      auto omega2_over_omega2_plus_rho = omega2 / (omega2 + rho);
      fm_eval_->eval(Gm, T * omega2_over_omega2_plus_rho, mmax);

      using std::sqrt;
      auto ooversqrto2prho_exp_2mplus1 = sqrt(omega2_over_omega2_plus_rho);
      for (auto m = 0; m <= mmax;
           ++m, ooversqrto2prho_exp_2mplus1 *= omega2_over_omega2_plus_rho) {
        Gm[m] *= ooversqrto2prho_exp_2mplus1;
      }
    } else {
      std::fill(Gm, Gm + mmax + 1, Real{0});
    }
  }

 private:
  std::shared_ptr<const FmEvalType> fm_eval_;  // need for odd K
};

/// core integral evaluator for \f$ ( \lambda \mathrm{erf}(\omega r) + \sigma
/// \mathrm{erfc}(\omega r) ) / r \f$ kernel
/// @note need extra scratch for Boys function values,
///       since need to call Boys engine twice
template <typename Real>
struct erfx_coulomb_gm_eval
    : private detail::CoreEvalScratch<erfx_coulomb_gm_eval<Real>> {
  typedef detail::CoreEvalScratch<erfx_coulomb_gm_eval<Real>> base_type;
  typedef Real value_type;

#ifndef LIBINT_USER_DEFINED_REAL
  using FmEvalType = libint2::FmEval_Chebyshev7<double>;
#else
  using FmEvalType = libint2::FmEval_Reference<scalar_type>;
#endif

  erfx_coulomb_gm_eval(unsigned int mmax, Real precision) : base_type(mmax) {
    fm_eval_ = FmEvalType::instance(mmax, precision);
  }
  void operator()(Real* Gm, Real rho, Real T, int mmax, Real omega, Real lambda,
                  Real sigma) {
    // \lambda \mathrm{erf}(\omega r) + \sigma \mathrm{erfc}(\omega r) = \sigma
    // - (\sigma - \lambda) \mathrm{erf}(\omega r)
    if (sigma != 0) {
      fm_eval_->eval(&base_type::Fm_[0], T, mmax);
      for (auto m = 0; m <= mmax; ++m) Gm[m] = sigma * base_type::Fm_[m];
    } else {
      std::fill(Gm, Gm + mmax + 1, Real{0});
    }
    const auto sigma_minus_lambda = sigma - lambda;
    if (omega > 0 && sigma_minus_lambda != 0) {
      auto omega2 = omega * omega;
      auto omega2_over_omega2_plus_rho = omega2 / (omega2 + rho);
      fm_eval_->eval(&base_type::Fm_[0], T * omega2_over_omega2_plus_rho, mmax);

      using std::sqrt;
      auto ooversqrto2prho_exp_2mplus1 =
          sqrt(omega2_over_omega2_plus_rho) * sigma_minus_lambda;
      for (auto m = 0; m <= mmax;
           ++m, ooversqrto2prho_exp_2mplus1 *= omega2_over_omega2_plus_rho) {
        Gm[m] -= ooversqrto2prho_exp_2mplus1 * base_type::Fm_[m];
      }
    }
  }

 private:
  std::shared_ptr<const FmEvalType> fm_eval_;  // need for odd K
};

}  // namespace os_core_ints

/*
 *  Slater geminal fitting is available only if have LAPACK
 */
#if HAVE_LAPACK
/*
f[x_] := - Exp[-\[Zeta] x] / \[Zeta];

ff[cc_, aa_, x_] := Sum[cc[[i]]*Exp[-aa[[i]] x^2], {i, 1, n}];
*/
template <typename Real>
Real fstg(Real zeta, Real x) {
  return -std::exp(-zeta * x) / zeta;
}

template <typename Real>
Real fngtg(const std::vector<Real>& cc, const std::vector<Real>& aa, Real x) {
  Real value = 0.0;
  const Real x2 = x * x;
  const unsigned int n = cc.size();
  for (unsigned int i = 0; i < n; ++i) value += cc[i] * std::exp(-aa[i] * x2);
  return value;
}

// --- weighting functions ---
// L2 error is weighted by ww(x)
// hence error is weighted by sqrt(ww(x))
template <typename Real>
Real wwtewklopper(Real x) {
  const Real x2 = x * x;
  return x2 * std::exp(-2 * x2);
}
template <typename Real>
Real wwcusp(Real x) {
  const Real x2 = x * x;
  const Real x6 = x2 * x2 * x2;
  return std::exp(-0.005 * x6);
}
// default is Tew-Klopper
template <typename Real>
Real ww(Real x) {
  // return wwtewklopper(x);
  return wwcusp(x);
}

template <typename Real>
Real norm(const std::vector<Real>& vec) {
  Real value = 0.0;
  const unsigned int n = vec.size();
  for (unsigned int i = 0; i < n; ++i) value += vec[i] * vec[i];
  return value;
}

template <typename Real>
void LinearSolveDamped(const std::vector<Real>& A, const std::vector<Real>& b,
                       Real lambda, std::vector<Real>& x) {
  const size_t n = b.size();
  std::vector<Real> Acopy(A);
  for (size_t m = 0; m < n; ++m) Acopy[m * n + m] *= (1 + lambda);
  std::vector<Real> e(b);

  // int info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, 1, &Acopy[0], n, &ipiv[0],
  // &e[0], n );
  {
    std::vector<int> ipiv(n);
    int n = b.size();
    int one = 1;
    int info;
    dgesv_(&n, &one, &Acopy[0], &n, &ipiv[0], &e[0], &n, &info);
    assert(info == 0);
  }

  x = e;
}

/**
 * computes a least-squares fit of \f$ -exp(-\zeta r_{12})/\zeta =
 * \sum_{i=1}^n c_i exp(-a_i r_{12}^2) \f$ on \f$ r_{12} \in [0, x_{\rm max}]
 * \f$ discretized to npts.
 * @param[in] n
 * @param[in] zeta
 * @param[out] geminal
 * @param[in] xmin
 * @param[in] xmax
 * @param[in] npts
 */
template <typename Real>
void stg_ng_fit(unsigned int n, Real zeta,
                std::vector<std::pair<Real, Real>>& geminal, Real xmin = 0.0,
                Real xmax = 10.0, unsigned int npts = 1001) {
  // initial guess
  std::vector<Real> cc(n, 1.0);  // coefficients
  std::vector<Real> aa(n);       // exponents
  for (unsigned int i = 0; i < n; ++i)
    aa[i] = std::pow(3.0, (i + 2 - (n + 1) / 2.0));

  // first rescale cc for ff[x] to match the norm of f[x]
  Real ffnormfac = 0.0;
  using std::sqrt;
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      ffnormfac += cc[i] * cc[j] / sqrt(aa[i] + aa[j]);
  const Real Nf = sqrt(2.0 * zeta) * zeta;
  const Real Nff = sqrt(2.0) / (sqrt(ffnormfac) * sqrt(sqrt(M_PI)));
  for (unsigned int i = 0; i < n; ++i) cc[i] *= -Nff / Nf;

  Real lambda0 = 1000;  // damping factor is initially set to 1000, eventually
                        // should end up at 0
  const Real nu = 3.0;  // increase/decrease the damping factor scale it by this
  const Real epsilon = 1e-15;  // convergence
  const unsigned int maxniter = 200;

  // grid points on which we will fit
  std::vector<Real> xi(npts);
  for (unsigned int i = 0; i < npts; ++i)
    xi[i] = xmin + (xmax - xmin) * i / (npts - 1);

  std::vector<Real> err(npts);

  const size_t nparams =
      2 * n;  // params = expansion coefficients + gaussian exponents
  std::vector<Real> J(npts * nparams);
  std::vector<Real> delta(nparams);

  //    std::cout << "iteration 0" << std::endl;
  //    for(unsigned int i=0; i<n; ++i)
  //      std::cout << cc[i] << " " << aa[i] << std::endl;

  Real errnormI;
  Real errnormIm1 = 1e3;
  bool converged = false;
  unsigned int iter = 0;
  while (!converged && iter < maxniter) {
    //      std::cout << "Iteration " << ++iter << ": lambda = " << lambda0/nu
    //      << std::endl;

    for (unsigned int i = 0; i < npts; ++i) {
      const Real x = xi[i];
      err[i] = (fstg(zeta, x) - fngtg(cc, aa, x)) * sqrt(ww(x));
    }
    errnormI = norm(err) / sqrt((Real)npts);

    //        std::cout << "|err|=" << errnormI << std::endl;
    converged = std::abs((errnormI - errnormIm1) / errnormIm1) <= epsilon;
    if (converged) break;
    errnormIm1 = errnormI;

    for (unsigned int i = 0; i < npts; ++i) {
      const Real x2 = xi[i] * xi[i];
      const Real sqrt_ww_x = sqrt(ww(xi[i]));
      const unsigned int ioffset = i * nparams;
      for (unsigned int j = 0; j < n; ++j)
        J[ioffset + j] = (std::exp(-aa[j] * x2)) * sqrt_ww_x;
      const unsigned int ioffsetn = ioffset + n;
      for (unsigned int j = 0; j < n; ++j)
        J[ioffsetn + j] = -sqrt_ww_x * x2 * cc[j] * std::exp(-aa[j] * x2);
    }

    std::vector<Real> A(nparams * nparams);
    for (size_t r = 0, rc = 0; r < nparams; ++r) {
      for (size_t c = 0; c < nparams; ++c, ++rc) {
        double Arc = 0.0;
        for (size_t i = 0, ir = r, ic = c; i < npts;
             ++i, ir += nparams, ic += nparams)
          Arc += J[ir] * J[ic];
        A[rc] = Arc;
      }
    }

    std::vector<Real> b(nparams);
    for (size_t r = 0; r < nparams; ++r) {
      Real br = 0.0;
      for (size_t i = 0, ir = r; i < npts; ++i, ir += nparams)
        br += J[ir] * err[i];
      b[r] = br;
    }

    // try decreasing damping first
    // if not successful try increasing damping until it results in a decrease
    // in the error
    lambda0 /= nu;
    for (int l = -1; l < 1000; ++l) {
      LinearSolveDamped(A, b, lambda0, delta);

      std::vector<double> cc_0(cc);
      for (unsigned int i = 0; i < n; ++i) cc_0[i] += delta[i];
      std::vector<double> aa_0(aa);
      for (unsigned int i = 0; i < n; ++i) aa_0[i] += delta[i + n];

      // if any of the exponents are negative the step is too large and need
      // to increase damping
      bool step_too_large = false;
      for (unsigned int i = 0; i < n; ++i)
        if (aa_0[i] < 0.0) {
          step_too_large = true;
          break;
        }
      if (!step_too_large) {
        std::vector<double> err_0(npts);
        for (unsigned int i = 0; i < npts; ++i) {
          const double x = xi[i];
          err_0[i] = (fstg(zeta, x) - fngtg(cc_0, aa_0, x)) * sqrt(ww(x));
        }
        const double errnorm_0 = norm(err_0) / sqrt((double)npts);
        if (errnorm_0 < errnormI) {
          cc = cc_0;
          aa = aa_0;
          break;
        } else  // step lead to increase of the error -- try dampening a bit
                // more
          lambda0 *= nu;
      } else  // too large of a step
        lambda0 *= nu;
    }  // done adjusting the damping factor

  }  // end of iterative minimization

  // if reached max # of iterations throw if the error is too terrible
  assert(not(iter == maxniter && errnormI > 1e-10));

  for (unsigned int i = 0; i < n; ++i)
    geminal[i] = std::make_pair(aa[i], cc[i]);
}
#endif

}  // end of namespace libint2

#endif  // C++ only
#endif  // header guard
