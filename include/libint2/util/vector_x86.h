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

#ifndef _libint2_src_lib_libint_vectorx86_h_
#define _libint2_src_lib_libint_vectorx86_h_

#include <libint2/util/cxxstd.h>
#include <libint2/util/type_traits.h>

#include <cmath>
#include <cstring>
#include <iostream>

#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__SSE2__) || defined(__SSE__) || defined(__AVX__)
#include <x86intrin.h>
#endif

#ifdef __SSE2__

namespace libint2 {
namespace simd {

/**
 * SIMD vector of 2 double-precision floating-point real numbers, operations on
 * which use SSE2 instructions available on all recent x86 hardware
 */
struct VectorSSEDouble {
  typedef double T;
  __m128d d;

  /**
   * creates a vector of default-initialized values.
   */
  VectorSSEDouble() {}

  /** Initializes all elements to the same value
   *  @param a the value to which all elements will be set
   */
  VectorSSEDouble(T a) { d = _mm_set_pd(a, a); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorSSEDouble(T (&a)[2]) { d = _mm_loadu_pd(&a[0]); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorSSEDouble(T a0, T a1) { d = _mm_set_pd(a1, a0); }

  /**
   * converts a 128-bit SSE double vector type to VectorSSEDouble
   */
  VectorSSEDouble(__m128d a) { d = a; }

  VectorSSEDouble& operator=(T a) {
    d = _mm_set_pd(a, a);
    return *this;
  }

  VectorSSEDouble& operator+=(VectorSSEDouble a) {
    d = _mm_add_pd(d, a.d);
    return *this;
  }

  VectorSSEDouble& operator-=(VectorSSEDouble a) {
    d = _mm_sub_pd(d, a.d);
    return *this;
  }

  VectorSSEDouble operator-() const {
    // TODO improve using bitflips
    const static __m128d minus_one = _mm_set_pd(-1.0, -1.0);
    ;
    VectorSSEDouble result;
    result.d = _mm_mul_pd(this->d, minus_one);
    return result;
  }

#if LIBINT2_CPLUSPLUS_STD >= 2011
  explicit
#endif
  operator double() const {
    double d0[2];
    ::memcpy(&(d0[0]), &d, sizeof(__m128d));
    // this segfaults on OS X
    //_mm_storeu_pd(&(d0[0]), d);
    //        // aligned alternative requires C++11's alignas, but efficiency
    //        should not matter here alignas(__m128d) double d0[2];
    //        _mm_store_pd(&(d0[0]), d);
    return d0[0];
  }

  /// implicit conversion to SSE 128-bit "register"
  operator __m128d() const { return d; }

  /// loads \c a to \c this
  void load(T const* a) { d = _mm_loadu_pd(a); }
  /// loads \c a to \c this  \sa load()
  /// @note \c a must be aligned to 16 bytes
  void load_aligned(T const* a) { d = _mm_load_pd(a); }
  /// writes \c this to \c a
  void convert(T* a) const { _mm_storeu_pd(&a[0], d); }
  /// writes \c this to \c a
  /// @note \c a must be aligned to 16 bytes
  void convert_aligned(T* a) const { _mm_store_pd(&a[0], d); }
};

//@{ arithmetic operators
inline VectorSSEDouble operator*(double a, VectorSSEDouble b) {
  VectorSSEDouble c;
  VectorSSEDouble _a(a);
  c.d = _mm_mul_pd(_a.d, b.d);
  return c;
}

inline VectorSSEDouble operator*(VectorSSEDouble a, double b) {
  VectorSSEDouble c;
  VectorSSEDouble _b(b);
  c.d = _mm_mul_pd(a.d, _b.d);
  return c;
}

inline VectorSSEDouble operator*(int a, VectorSSEDouble b) {
  if (a == 1)
    return b;
  else {
    VectorSSEDouble c;
    VectorSSEDouble _a((double)a);
    c.d = _mm_mul_pd(_a.d, b.d);
    return c;
  }
}

inline VectorSSEDouble operator*(VectorSSEDouble a, int b) {
  if (b == 1)
    return a;
  else {
    VectorSSEDouble c;
    VectorSSEDouble _b((double)b);
    c.d = _mm_mul_pd(a.d, _b.d);
    return c;
  }
}

inline VectorSSEDouble operator*(VectorSSEDouble a, VectorSSEDouble b) {
  VectorSSEDouble c;
  c.d = _mm_mul_pd(a.d, b.d);
  return c;
}

inline VectorSSEDouble operator+(VectorSSEDouble a, VectorSSEDouble b) {
  VectorSSEDouble c;
  c.d = _mm_add_pd(a.d, b.d);
  return c;
}

inline VectorSSEDouble operator-(VectorSSEDouble a, VectorSSEDouble b) {
  VectorSSEDouble c;
  c.d = _mm_sub_pd(a.d, b.d);
  return c;
}

inline VectorSSEDouble operator/(VectorSSEDouble a, VectorSSEDouble b) {
  VectorSSEDouble c;
  c.d = _mm_div_pd(a.d, b.d);
  return c;
}

#if defined(__FMA__)
inline VectorSSEDouble fma_plus(VectorSSEDouble a, VectorSSEDouble b,
                                VectorSSEDouble c) {
  VectorSSEDouble d;
  d.d = _mm_fmadd_pd(a.d, b.d, c.d);
  return d;
}
inline VectorSSEDouble fma_minus(VectorSSEDouble a, VectorSSEDouble b,
                                 VectorSSEDouble c) {
  VectorSSEDouble d;
  d.d = _mm_fmsub_pd(a.d, b.d, c.d);
  return d;
}
#elif defined(__FMA4__)
inline VectorSSEDouble fma_plus(VectorSSEDouble a, VectorSSEDouble b,
                                VectorSSEDouble c) {
  VectorSSEDouble d;
  d.d = _mm_macc_pd(a.d, b.d, c.d);
  return d;
}
inline VectorSSEDouble fma_minus(VectorSSEDouble a, VectorSSEDouble b,
                                 VectorSSEDouble c) {
  VectorSSEDouble d;
  d.d = _mm_msub_pd(a.d, b.d, c.d);
  return d;
}
#endif

/// Horizontal add
/// @param a input vector = {a[0], a[1]}
/// @return a[0] + a[1]
inline double horizontal_add(VectorSSEDouble const& a) {
//      Agner Fog's version
#if defined(__SSE3__)
  __m128d t1 = _mm_hadd_pd(a, a);
  return _mm_cvtsd_f64(t1);
#else  // SSE2 only
  __m128 t0 = _mm_castpd_ps(a);
  __m128d t1 = _mm_castps_pd(_mm_movehl_ps(t0, t0));
  __m128d t2 = _mm_add_sd(a, t1);
  return _mm_cvtsd_f64(t2);
#endif
}

/// Horizontal add of a pair of vectors
/// @param a input vector = {a[0], a[1]}
/// @param b input vector = {b[0], b[1]}
/// @return {a[0] + a[1], b[0] + b[1]}
inline VectorSSEDouble horizontal_add(VectorSSEDouble const& a,
                                      VectorSSEDouble const& b) {
#if defined(__SSE3__)
  return _mm_hadd_pd(a, b);
#else  // will be very inefficient without SSE3
  return VectorSSEDouble(horizontal_add(a), horizontal_add(b));
#endif
}

//@}

//@{ standard functions
inline VectorSSEDouble exp(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
  VectorSSEDouble result;
  result.d = _mm_exp_pd(a.d);
#else
  double a_d[2];
  a.convert(a_d);
  for (int i = 0; i < 2; ++i) a_d[i] = std::exp(a_d[i]);
  VectorSSEDouble result(a_d);
#endif
  return result;
}
inline VectorSSEDouble sqrt(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
  VectorSSEDouble result;
  result.d = _mm_sqrt_pd(a.d);
#else
  double a_d[2];
  a.convert(a_d);
  for (int i = 0; i < 2; ++i) a_d[i] = std::sqrt(a_d[i]);
  VectorSSEDouble result(a_d);
#endif
  return result;
}
inline VectorSSEDouble erf(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
  VectorSSEDouble result;
  result.d = _mm_erf_pd(a.d);
#else
  double a_d[2];
  a.convert(a_d);
  for (int i = 0; i < 2; ++i) a_d[i] = ::erf(a_d[i]);
  VectorSSEDouble result(a_d);
#endif
  return result;
}
inline VectorSSEDouble erfc(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
  VectorSSEDouble result;
  result.d = _mm_erfc_pd(a.d);
#else
  double a_d[2];
  a.convert(a_d);
  for (int i = 0; i < 2; ++i) a_d[i] = ::erfc(a_d[i]);
  VectorSSEDouble result(a_d);
#endif
  return result;
}
//@}

};  // namespace simd
};  // namespace libint2

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os,
                                libint2::simd::VectorSSEDouble a) {
  double ad[2];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "}";
  return os;
}
//@}

namespace libint2 {

//@{ vector traits of VectorSSEDouble

template <>
struct is_vector<simd::VectorSSEDouble> {
  static const bool value = true;
};

template <>
struct vector_traits<simd::VectorSSEDouble> {
  typedef double scalar_type;
  static const size_t extent = 2;
};

//@}

}  // namespace libint2

#endif  // SSE2-only

#ifdef __SSE__

namespace libint2 {
namespace simd {

/**
 * SIMD vector of 4 single-precision floating-point real numbers, operations on
 * which use SSE instructions available on all recent x86 hardware.
 */
struct VectorSSEFloat {
  typedef float T;
  __m128 d;

  /**
   * creates a vector of default-initialized values.
   */
  VectorSSEFloat() {}

  /** Initializes all elements to the same value
   *  @param a the value to which all elements will be set
   */
  VectorSSEFloat(T a) { d = _mm_set_ps(a, a, a, a); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorSSEFloat(T (&a)[4]) { d = _mm_loadu_ps(&a[0]); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorSSEFloat(T a0, T a1, T a2, T a3) { d = _mm_set_ps(a3, a2, a1, a0); }

  VectorSSEFloat& operator=(T a) {
    d = _mm_set_ps(a, a, a, a);
    return *this;
  }

  VectorSSEFloat& operator+=(VectorSSEFloat a) {
    d = _mm_add_ps(d, a.d);
    return *this;
  }

  VectorSSEFloat& operator-=(VectorSSEFloat a) {
    d = _mm_sub_ps(d, a.d);
    return *this;
  }

  VectorSSEFloat operator-() const {
    // TODO improve using bitflips
    const static __m128 minus_one = _mm_set_ps(-1.0, -1.0, -1.0, -1.0);
    ;
    VectorSSEFloat result;
    result.d = _mm_mul_ps(this->d, minus_one);
    return result;
  }

#if LIBINT2_CPLUSPLUS_STD >= 2011
  explicit
#endif
  operator float() const {
    float d0[4];
    ::memcpy(&(d0[0]), &d, sizeof(__m128));
    // this segfaults on OS X
    //_mm_storeu_ps(&(d0[0]), d);
    //        // aligned alternative requires C++11's alignas, but efficiency
    //        should not matter here alignas(__m128) float d0[4];
    //        _mm_store_ps(&(d0[0]), d);
    return d0[0];
  }

#if LIBINT2_CPLUSPLUS_STD >= 2011
  explicit
#endif
  operator double() const {
    const float result_flt = this->operator float();
    return result_flt;
  }

  /// implicit conversion to SSE 128-bit "register"
  operator __m128() const { return d; }

  /// loads \c a to \c this
  void load(T const* a) { d = _mm_loadu_ps(a); }
  /// loads \c a to \c this  \sa load()
  /// @note \c a must be aligned to 16 bytes
  void load_aligned(T const* a) { d = _mm_load_ps(a); }
  /// writes \c this to \c a
  void convert(T* a) const { _mm_storeu_ps(&a[0], d); }
  /// writes \c this to \c a
  /// @note \c a must be aligned to 32 bytes
  void convert_aligned(T* a) const { _mm_store_ps(&a[0], d); }
};

//@{ arithmetic operators
inline VectorSSEFloat operator*(float a, VectorSSEFloat b) {
  VectorSSEFloat c;
  VectorSSEFloat _a(a);
  c.d = _mm_mul_ps(_a.d, b.d);
  return c;
}

inline VectorSSEFloat operator*(VectorSSEFloat a, float b) {
  VectorSSEFloat c;
  VectorSSEFloat _b(b);
  c.d = _mm_mul_ps(a.d, _b.d);
  return c;
}

// narrows a!
inline VectorSSEFloat operator*(double a, VectorSSEFloat b) {
  VectorSSEFloat c;
  VectorSSEFloat _a((float)a);
  c.d = _mm_mul_ps(_a.d, b.d);
  return c;
}

// narrows b!
inline VectorSSEFloat operator*(VectorSSEFloat a, double b) {
  VectorSSEFloat c;
  VectorSSEFloat _b((float)b);
  c.d = _mm_mul_ps(a.d, _b.d);
  return c;
}

inline VectorSSEFloat operator*(int a, VectorSSEFloat b) {
  if (a == 1)
    return b;
  else {
    VectorSSEFloat c;
    VectorSSEFloat _a((float)a);
    c.d = _mm_mul_ps(_a.d, b.d);
    return c;
  }
}

inline VectorSSEFloat operator*(VectorSSEFloat a, int b) {
  if (b == 1)
    return a;
  else {
    VectorSSEFloat c;
    VectorSSEFloat _b((float)b);
    c.d = _mm_mul_ps(a.d, _b.d);
    return c;
  }
}

inline VectorSSEFloat operator*(VectorSSEFloat a, VectorSSEFloat b) {
  VectorSSEFloat c;
  c.d = _mm_mul_ps(a.d, b.d);
  return c;
}

inline VectorSSEFloat operator+(VectorSSEFloat a, VectorSSEFloat b) {
  VectorSSEFloat c;
  c.d = _mm_add_ps(a.d, b.d);
  return c;
}

inline VectorSSEFloat operator-(VectorSSEFloat a, VectorSSEFloat b) {
  VectorSSEFloat c;
  c.d = _mm_sub_ps(a.d, b.d);
  return c;
}

inline VectorSSEFloat operator/(VectorSSEFloat a, VectorSSEFloat b) {
  VectorSSEFloat c;
  c.d = _mm_div_ps(a.d, b.d);
  return c;
}

#if defined(__FMA__)
inline VectorSSEFloat fma_plus(VectorSSEFloat a, VectorSSEFloat b,
                               VectorSSEFloat c) {
  VectorSSEFloat d;
  d.d = _mm_fmadd_ps(a.d, b.d, c.d);
  return d;
}
inline VectorSSEFloat fma_minus(VectorSSEFloat a, VectorSSEFloat b,
                                VectorSSEFloat c) {
  VectorSSEFloat d;
  d.d = _mm_fmsub_ps(a.d, b.d, c.d);
  return d;
}
#elif defined(__FMA4__)
inline VectorSSEFloat fma_plus(VectorSSEFloat a, VectorSSEFloat b,
                               VectorSSEFloat c) {
  VectorSSEFloat d;
  d.d = _mm_macc_ps(a.d, b.d, c.d);
  return d;
}
inline VectorSSEFloat fma_minus(VectorSSEFloat a, VectorSSEFloat b,
                                VectorSSEFloat c) {
  VectorSSEFloat d;
  d.d = _mm_msub_ps(a.d, b.d, c.d);
  return d;
}
#endif

//@}

//@{ standard functions
inline VectorSSEFloat exp(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
  VectorSSEFloat result;
  result.d = _mm_exp_ps(a.d);
#else
  float a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = std::exp(a_d[i]);
  VectorSSEFloat result(a_d);
#endif
  return result;
}
inline VectorSSEFloat sqrt(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
  VectorSSEFloat result;
  result.d = _mm_sqrt_ps(a.d);
#else
  float a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = std::sqrt(a_d[i]);
  VectorSSEFloat result(a_d);
#endif
  return result;
}
inline VectorSSEFloat erf(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
  VectorSSEFloat result;
  result.d = _mm_erf_ps(a.d);
#else
  float a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::erf(a_d[i]);
  VectorSSEFloat result(a_d);
#endif
  return result;
}
inline VectorSSEFloat erfc(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
  VectorSSEFloat result;
  result.d = _mm_erfc_ps(a.d);
#else
  float a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::erfc(a_d[i]);
  VectorSSEFloat result(a_d);
#endif
  return result;
}
//@}

};  // namespace simd
};  // namespace libint2

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os,
                                libint2::simd::VectorSSEFloat a) {
  float ad[4];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "," << ad[2] << "," << ad[3] << "}";
  return os;
}
//@}

namespace libint2 {

//@{ vector traits of VectorSSEFloat

template <>
struct is_vector<simd::VectorSSEFloat> {
  static const bool value = true;
};

template <>
struct vector_traits<simd::VectorSSEFloat> {
  typedef float scalar_type;
  static const size_t extent = 4;
};

//@}

}  // namespace libint2

#endif  // SSE-only

#ifdef __AVX__

namespace libint2 {
namespace simd {

/**
 * SIMD vector of 4 double-precision floating-point real numbers, operations on
 * which use AVX instructions available on recent x86 hardware from Intel
 * (starting with Sandy Bridge processors released in 2011) and AMD (starting
 * with Bulldozer released in 2011).
 */
struct VectorAVXDouble {
  typedef double T;
  __m256d d;

  /**
   * creates a vector of default-initialized values.
   */
  VectorAVXDouble() {}

  /** Initializes all elements to the same value
   *  @param a the value to which all elements will be set
   */
  VectorAVXDouble(T a) { d = _mm256_set_pd(a, a, a, a); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorAVXDouble(T (&a)[4]) { d = _mm256_loadu_pd(&a[0]); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorAVXDouble(T a0, T a1, T a2, T a3) { d = _mm256_set_pd(a3, a2, a1, a0); }

  /**
   * converts a 256-bit AVX double vector type to VectorAVXDouble
   */
  VectorAVXDouble(__m256d a) { d = a; }

  VectorAVXDouble& operator=(T a) {
    d = _mm256_set_pd(a, a, a, a);
    return *this;
  }

  VectorAVXDouble& operator+=(VectorAVXDouble a) {
    d = _mm256_add_pd(d, a.d);
    return *this;
  }

  VectorAVXDouble& operator-=(VectorAVXDouble a) {
    d = _mm256_sub_pd(d, a.d);
    return *this;
  }

  VectorAVXDouble operator-() const {
    // TODO improve using bitflips
    const static __m256d minus_one = _mm256_set_pd(-1.0, -1.0, -1.0, -1.0);
    ;
    VectorAVXDouble result;
    result.d = _mm256_mul_pd(this->d, minus_one);
    return result;
  }

#if LIBINT2_CPLUSPLUS_STD >= 2011
  explicit
#endif
  operator double() const {
    double d0[4];
    ::memcpy(&(d0[0]), &d, sizeof(__m256d));
    // this segfaults on OS X
    //        _mm256_storeu_pd(&d0[0], d);
    //        // aligned alternative requires C++11's alignas, but efficiency
    //        should not matter here
    //        //        alignas(__m256d) double d0[4];
    //        //        _mm256_store_pd(&(d0[0]), d);
    return d0[0];
  }

  /// implicit conversion to AVX 256-bit "register"
  operator __m256d() const { return d; }

  /// loads \c a to \c this
  void load(T const* a) { d = _mm256_loadu_pd(a); }
  /// loads \c a to \c this  \sa load()
  /// @note \c a must be aligned to 32 bytes
  void load_aligned(T const* a) { d = _mm256_load_pd(a); }
  /// writes \c this to \c a
  void convert(T* a) const { _mm256_storeu_pd(&a[0], d); }
  /// writes \c this to \c a
  /// @note \c a must be aligned to 32 bytes
  void convert_aligned(T* a) const { _mm256_store_pd(&a[0], d); }
};

//@{ arithmetic operators
inline VectorAVXDouble operator*(double a, VectorAVXDouble b) {
  VectorAVXDouble c;
  VectorAVXDouble _a(a);
  c.d = _mm256_mul_pd(_a.d, b.d);
  return c;
}

inline VectorAVXDouble operator*(VectorAVXDouble a, double b) {
  VectorAVXDouble c;
  VectorAVXDouble _b(b);
  c.d = _mm256_mul_pd(a.d, _b.d);
  return c;
}

inline VectorAVXDouble operator*(int a, VectorAVXDouble b) {
  if (a == 1)
    return b;
  else {
    VectorAVXDouble c;
    VectorAVXDouble _a((double)a);
    c.d = _mm256_mul_pd(_a.d, b.d);
    return c;
  }
}

inline VectorAVXDouble operator*(VectorAVXDouble a, int b) {
  if (b == 1)
    return a;
  else {
    VectorAVXDouble c;
    VectorAVXDouble _b((double)b);
    c.d = _mm256_mul_pd(a.d, _b.d);
    return c;
  }
}

inline VectorAVXDouble operator*(VectorAVXDouble a, VectorAVXDouble b) {
  VectorAVXDouble c;
  c.d = _mm256_mul_pd(a.d, b.d);
  return c;
}

inline VectorAVXDouble operator+(VectorAVXDouble a, VectorAVXDouble b) {
  VectorAVXDouble c;
  c.d = _mm256_add_pd(a.d, b.d);
  return c;
}

inline VectorAVXDouble operator+(int a, VectorAVXDouble b) {
  if (a == 0) {
    return b;
  } else {
    VectorAVXDouble c;
    VectorAVXDouble _a = (static_cast<double>(a));
    c.d = _mm256_add_pd(_a.d, b.d);
    return c;
  }
}

inline VectorAVXDouble operator+(VectorAVXDouble a, int b) {
  if (b == 0) {
    return a;
  } else {
    VectorAVXDouble c;
    VectorAVXDouble _b = (static_cast<double>(b));
    c.d = _mm256_add_pd(a.d, _b.d);
    return c;
  }
}

inline VectorAVXDouble operator-(VectorAVXDouble a, VectorAVXDouble b) {
  VectorAVXDouble c;
  c.d = _mm256_sub_pd(a.d, b.d);
  return c;
}

inline VectorAVXDouble operator/(VectorAVXDouble a, VectorAVXDouble b) {
  VectorAVXDouble c;
  c.d = _mm256_div_pd(a.d, b.d);
  return c;
}

#if defined(__FMA__)
inline VectorAVXDouble fma_plus(VectorAVXDouble a, VectorAVXDouble b,
                                VectorAVXDouble c) {
  VectorAVXDouble d;
  d.d = _mm256_fmadd_pd(a.d, b.d, c.d);
  return d;
}
inline VectorAVXDouble fma_minus(VectorAVXDouble a, VectorAVXDouble b,
                                 VectorAVXDouble c) {
  VectorAVXDouble d;
  d.d = _mm256_fmsub_pd(a.d, b.d, c.d);
  return d;
}
#elif defined(__FMA4__)
inline VectorAVXDouble fma_plus(VectorAVXDouble a, VectorAVXDouble b,
                                VectorAVXDouble c) {
  VectorAVXDouble d;
  d.d = _mm256_macc_pd(a.d, b.d, c.d);
  return d;
}
inline VectorAVXDouble fma_minus(VectorAVXDouble a, VectorAVXDouble b,
                                 VectorAVXDouble c) {
  VectorAVXDouble d;
  d.d = _mm256_msub_pd(a.d, b.d, c.d);
  return d;
}
#endif

/// Horizontal add
/// @param a input vector = {a[0], a[1], a[2], a[3]}
/// @return a[0] + a[1] + a[2] + a[3]
inline double horizontal_add(VectorAVXDouble const& a) {
  __m256d s = _mm256_hadd_pd(a, a);
  return ((double*)&s)[0] + ((double*)&s)[2];
  //      Agner Fog's version
  //      __m256d t1 = _mm256_hadd_pd(a,a);
  //      __m128d t2 = _mm256_extractf128_pd(t1,1);
  //      __m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1),t2);
  //      return _mm_cvtsd_f64(t3);
}

/// Horizontal add of a pair of vectors
/// @param a input vector = {a[0], a[1], a[2], a[3]}
/// @param b input vector = {b[0], b[1], b[2], b[3]}
/// @return {a[0] + a[1] + a[2] + a[3], b[0] + b[1] + b[2] + b[3]}
inline VectorSSEDouble horizontal_add(VectorAVXDouble const& a,
                                      VectorAVXDouble const& b) {
  // solution from
  // http://stackoverflow.com/questions/9775538/fastest-way-to-do-horizontal-vector-sum-with-avx-instructions
  __m256d sum = _mm256_hadd_pd(a, b);
  // extract upper 128 bits of result
  __m128d sum_high = _mm256_extractf128_pd(sum, 1);
  // add upper 128 bits of sum to its lower 128 bits
  return _mm_add_pd(sum_high, _mm256_castpd256_pd128(sum));
}

/// Horizontal add of a set of 4 vectors
/// @param a input vector = {a[0], a[1], a[2], a[3]}
/// @param b input vector = {b[0], b[1], b[2], b[3]}
/// @param c input vector = {c[0], c[1], c[2], c[3]}
/// @param d input vector = {d[0], d[1], d[2], d[3]}
/// @return {a[0] + a[1] + a[2] + a[3], b[0] + b[1] + b[2] + b[3], c[0] + c[1] +
/// c[2] + c[3], d[0] + d[1] + d[2] + d[3]}
inline VectorAVXDouble horizontal_add(VectorAVXDouble const& a,
                                      VectorAVXDouble const& b,
                                      VectorAVXDouble const& c,
                                      VectorAVXDouble const& d) {
  // solution from
  // http://stackoverflow.com/questions/10833234/4-horizontal-double-precision-sums-in-one-go-with-avx?lq=1
  // {a[0]+a[1], b[0]+b[1], a[2]+a[3], b[2]+b[3]}
  __m256d sumab = _mm256_hadd_pd(a, b);
  // {c[0]+c[1], d[0]+d[1], c[2]+c[3], d[2]+d[3]}
  __m256d sumcd = _mm256_hadd_pd(c, d);

  // {a[0]+a[1], b[0]+b[1], c[2]+c[3], d[2]+d[3]}
  __m256d blend = _mm256_blend_pd(sumab, sumcd, 0b1100);
  // {a[2]+a[3], b[2]+b[3], c[0]+c[1], d[0]+d[1]}
  __m256d perm = _mm256_permute2f128_pd(sumab, sumcd, 0x21);

  return _mm256_add_pd(perm, blend);
}

//@}

//@{ standard functions
inline VectorAVXDouble exp(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
  VectorAVXDouble result;
  result.d = _mm256_exp_pd(a.d);
#else
  double a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::exp(a_d[i]);
  VectorAVXDouble result(a_d);
#endif
  return result;
}
inline VectorAVXDouble sqrt(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
  VectorAVXDouble result;
  result.d = _mm256_sqrt_pd(a.d);
#else
  double a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::sqrt(a_d[i]);
  VectorAVXDouble result(a_d);
#endif
  return result;
}
inline VectorAVXDouble erf(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
  VectorAVXDouble result;
  result.d = _mm256_erf_pd(a.d);
#else
  double a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::erf(a_d[i]);
  VectorAVXDouble result(a_d);
#endif
  return result;
}
inline VectorAVXDouble erfc(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
  VectorAVXDouble result;
  result.d = _mm256_erfc_pd(a.d);
#else
  double a_d[4];
  a.convert(a_d);
  for (int i = 0; i < 4; ++i) a_d[i] = ::erfc(a_d[i]);
  VectorAVXDouble result(a_d);
#endif
  return result;
}

};  // namespace simd
};  // namespace libint2

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os,
                                libint2::simd::VectorAVXDouble a) {
  double ad[4];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "," << ad[2] << "," << ad[3] << "}";
  return os;
}

namespace libint2 {

//@{ vector traits of VectorAVXDouble

template <>
struct is_vector<simd::VectorAVXDouble> {
  static const bool value = true;
};

template <>
struct vector_traits<simd::VectorAVXDouble> {
  typedef double scalar_type;
  static const size_t extent = 4;
};

//@}

}  // namespace libint2

namespace libint2 {
namespace simd {

/**
 * SIMD vector of 8 single-precision floating-point real numbers, operations on
 * which use AVX instructions available on recent x86 hardware from Intel
 * (starting with Sandy Bridge processors released in 2011) and AMD (starting
 * with Bulldozer released in 2011).
 */
struct VectorAVXFloat {
  typedef float T;
  __m256 d;

  /**
   * creates a vector of default-initialized values.
   */
  VectorAVXFloat() {}

  /** Initializes all elements to the same value
   *  @param a the value to which all elements will be set
   */
  VectorAVXFloat(T a) { d = _mm256_set_ps(a, a, a, a, a, a, a, a); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorAVXFloat(T (&a)[8]) { d = _mm256_loadu_ps(&a[0]); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  VectorAVXFloat(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7) {
    d = _mm256_set_ps(a0, a1, a2, a3, a4, a5, a6, a7);
  }

  VectorAVXFloat& operator=(T a) {
    d = _mm256_set_ps(a, a, a, a, a, a, a, a);
    return *this;
  }

  VectorAVXFloat& operator+=(VectorAVXFloat a) {
    d = _mm256_add_ps(d, a.d);
    return *this;
  }

  VectorAVXFloat& operator-=(VectorAVXFloat a) {
    d = _mm256_sub_ps(d, a.d);
    return *this;
  }

  VectorAVXFloat operator-() const {
    // TODO improve using bitflips
    const static __m256 minus_one =
        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
    ;
    VectorAVXFloat result;
    result.d = _mm256_mul_ps(this->d, minus_one);
    return result;
  }

  explicit operator float() const {
    float d0[8];
    ::memcpy(&(d0[0]), &d, sizeof(__m256));
    // this segfaults on OS X
    //        _mm256_storeu_ps(&d0[0], d);
    //        // aligned alternative requires C++11's alignas, but efficiency
    //        should not matter here
    //        //        alignas(__m256) float d0[8];
    //        //        _mm256_store_ps(&(d0[0]), d);
    return d0[0];
  }

  explicit operator double() const {
    const float result_flt = this->operator float();
    return result_flt;
  }

  void convert(T (&a)[8]) const { _mm256_storeu_ps(&a[0], d); }
};

//@{ arithmetic operators
inline VectorAVXFloat operator*(double a, VectorAVXFloat b) {
  VectorAVXFloat c;
  VectorAVXFloat _a(a);
  c.d = _mm256_mul_ps(_a.d, b.d);
  return c;
}

inline VectorAVXFloat operator*(VectorAVXFloat a, double b) {
  VectorAVXFloat c;
  VectorAVXFloat _b(b);
  c.d = _mm256_mul_ps(a.d, _b.d);
  return c;
}

inline VectorAVXFloat operator*(int a, VectorAVXFloat b) {
  if (a == 1)
    return b;
  else {
    VectorAVXFloat c;
    VectorAVXFloat _a((float)a);
    c.d = _mm256_mul_ps(_a.d, b.d);
    return c;
  }
}

inline VectorAVXFloat operator*(VectorAVXFloat a, int b) {
  if (b == 1)
    return a;
  else {
    VectorAVXFloat c;
    VectorAVXFloat _b((float)b);
    c.d = _mm256_mul_ps(a.d, _b.d);
    return c;
  }
}

inline VectorAVXFloat operator*(VectorAVXFloat a, VectorAVXFloat b) {
  VectorAVXFloat c;
  c.d = _mm256_mul_ps(a.d, b.d);
  return c;
}

inline VectorAVXFloat operator+(VectorAVXFloat a, VectorAVXFloat b) {
  VectorAVXFloat c;
  c.d = _mm256_add_ps(a.d, b.d);
  return c;
}

inline VectorAVXFloat operator-(VectorAVXFloat a, VectorAVXFloat b) {
  VectorAVXFloat c;
  c.d = _mm256_sub_ps(a.d, b.d);
  return c;
}

inline VectorAVXFloat operator/(VectorAVXFloat a, VectorAVXFloat b) {
  VectorAVXFloat c;
  c.d = _mm256_div_ps(a.d, b.d);
  return c;
}

#if defined(__FMA__)
inline VectorAVXFloat fma_plus(VectorAVXFloat a, VectorAVXFloat b,
                               VectorAVXFloat c) {
  VectorAVXFloat d;
  d.d = _mm256_fmadd_ps(a.d, b.d, c.d);
  return d;
}
inline VectorAVXFloat fma_minus(VectorAVXFloat a, VectorAVXFloat b,
                                VectorAVXFloat c) {
  VectorAVXFloat d;
  d.d = _mm256_fmsub_ps(a.d, b.d, c.d);
  return d;
}
#elif defined(__FMA4__)
inline VectorAVXFloat fma_plus(VectorAVXFloat a, VectorAVXFloat b,
                               VectorAVXFloat c) {
  VectorAVXFloat d;
  d.d = _mm256_facc_ps(a.d, b.d, c.d);
  return d;
}
inline VectorAVXFloat fma_minus(VectorAVXFloat a, VectorAVXFloat b,
                                VectorAVXFloat c) {
  VectorAVXFloat d;
  d.d = _mm256_fsub_ps(a.d, b.d, c.d);
  return d;
}
#endif

//@}

//@{ standard functions
inline VectorAVXFloat exp(VectorAVXFloat a) {
#if HAVE_INTEL_SVML
  VectorAVXFloat result;
  result.d = _mm256_exp_ps(a.d);
#else
  float a_d[8];
  a.convert(a_d);
  for (int i = 0; i < 8; ++i) a_d[i] = ::exp(a_d[i]);
  VectorAVXFloat result(a_d);
#endif
  return result;
}
inline VectorAVXFloat sqrt(VectorAVXFloat a) {
#if HAVE_INTEL_SVML
  VectorAVXFloat result;
  result.d = _mm256_sqrt_ps(a.d);
#else
  float a_d[8];
  a.convert(a_d);
  for (int i = 0; i < 8; ++i) a_d[i] = ::sqrt(a_d[i]);
  VectorAVXFloat result(a_d);
#endif
  return result;
}
inline VectorAVXFloat erf(VectorAVXFloat a) {
#if HAVE_INTEL_SVML
  VectorAVXFloat result;
  result.d = _mm256_erf_ps(a.d);
#else
  float a_d[8];
  a.convert(a_d);
  for (int i = 0; i < 8; ++i) a_d[i] = ::erf(a_d[i]);
  VectorAVXFloat result(a_d);
#endif
  return result;
}
inline VectorAVXFloat erfc(VectorAVXFloat a) {
#if HAVE_INTEL_SVML
  VectorAVXFloat result;
  result.d = _mm256_erfc_ps(a.d);
#else
  float a_d[8];
  a.convert(a_d);
  for (int i = 0; i < 8; ++i) a_d[i] = ::erfc(a_d[i]);
  VectorAVXFloat result(a_d);
#endif
  return result;
}
//@}

};  // namespace simd
};  // namespace libint2

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os,
                                libint2::simd::VectorAVXFloat a) {
  float ad[8];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "," << ad[2] << "," << ad[3] << ","
     << ad[4] << "," << ad[5] << "," << ad[6] << "," << ad[7] << "}";
  return os;
}
//@}

namespace libint2 {

//@{ vector traits of VectorAVXFloat

template <>
struct is_vector<simd::VectorAVXFloat> {
  static const bool value = true;
};

template <>
struct vector_traits<simd::VectorAVXFloat> {
  typedef float scalar_type;
  static const size_t extent = 8;
};

//@}

}  // namespace libint2

#endif  // AVX-only

#ifdef LIBINT2_HAVE_AGNER_VECTORCLASS
#include <vectorclass.h>
#endif

#endif  // header guard
