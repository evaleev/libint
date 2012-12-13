
#ifndef _libint2_src_lib_libint_vectorx86_h_
#define _libint2_src_lib_libint_vectorx86_h_

#ifdef __SSE2__

#include <emmintrin.h>
#include <cmath>
#include <iostream>

namespace libint2 { namespace simd {

  /**
   * SIMD vector of 2 double-precision floating-point real numbers, operations on which use SSE2 instructions
   * available on all recent x86 hardware
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
      VectorSSEDouble(T a) {
        d = _mm_set_pd(a, a);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorSSEDouble(T (&a)[2]) {
        d = _mm_loadu_pd(&a[0]);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorSSEDouble(T a0, T a1) {
        d = _mm_set_pd(a0, a1);
      }

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
        const static __m128d minus_one = _mm_set_pd(-1.0, -1.0);;
        VectorSSEDouble result;
        result.d = _mm_mul_pd(this->d, minus_one);
        return result;
      }

//      operator double() const {
//        double d0[2];
//        _mm_store_sd(&(d0[0]), d);
//        return d0[0];
//      }

      void convert(double(&a)[2]) const {
        _mm_storeu_pd(&a[0], d);
      }

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
  inline VectorSSEDouble fma_plus(VectorSSEDouble a, VectorSSEDouble b, VectorSSEDouble c) {
    VectorSSEDouble d;
    d.d = _mm_fmadd_pd(a.d, b.d, c.d);
    return d;
  }
  inline VectorSSEDouble fma_minus(VectorSSEDouble a, VectorSSEDouble b, VectorSSEDouble c) {
    VectorSSEDouble d;
    d.d = _mm_fmsub_pd(a.d, b.d, c.d);
    return d;
  }
#elif defined(__FMA4__)
  inline VectorSSEDouble fma_plus(VectorSSEDouble a, VectorSSEDouble b, VectorSSEDouble c) {
    VectorSSEDouble d;
    d.d = _mm_macc_pd(a.d, b.d, c.d);
    return d;
  }
  inline VectorSSEDouble fma_minus(VectorSSEDouble a, VectorSSEDouble b, VectorSSEDouble c) {
    VectorSSEDouble d;
    d.d = _mm_msub_pd(a.d, b.d, c.d);
    return d;
  }
#endif

  //@}

  //@{ standard functions
  inline VectorSSEDouble exp(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
    VectorSSEDouble result;
    result.d = _mm_exp_pd(a.d);
#else
    double a_d[2]; a.convert(a_d);
    for(int i=0; i<2; ++i) a_d[i] = std::exp(a_d[i]);
    VectorSSEDouble result(a_d);
#endif
    return result;
  }
  inline VectorSSEDouble sqrt(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
    VectorSSEDouble result;
    result.d = _mm_sqrt_pd(a.d);
#else
    double a_d[2]; a.convert(a_d);
    for(int i=0; i<2; ++i) a_d[i] = std::sqrt(a_d[i]);
    VectorSSEDouble result(a_d);
#endif
    return result;
  }
  inline VectorSSEDouble erf(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
    VectorSSEDouble result;
    result.d = _mm_erf_pd(a.d);
#else
    double a_d[2]; a.convert(a_d);
    for(int i=0; i<2; ++i) a_d[i] = ::erf(a_d[i]);
    VectorSSEDouble result(a_d);
#endif
    return result;
  }
  inline VectorSSEDouble erfc(VectorSSEDouble a) {
#if HAVE_INTEL_SVML
    VectorSSEDouble result;
    result.d = _mm_erfc_pd(a.d);
#else
    double a_d[2]; a.convert(a_d);
    for(int i=0; i<2; ++i) a_d[i] = ::erfc(a_d[i]);
    VectorSSEDouble result(a_d);
#endif
    return result;
  }
  //@}

};}; // namespace libint2::simd

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os, libint2::simd::VectorSSEDouble a) {
  double ad[2];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "}";
  return os;
}
//@}

#endif // SSE2-only

#ifdef __SSE__

#include <xmmintrin.h>

namespace libint2 { namespace simd {

  /**
   * SIMD vector of 4 single-precision floating-point real numbers, operations on which use SSE instructions
   * available on all recent x86 hardware.
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
      VectorSSEFloat(T a) {
        d = _mm_set_ps(a, a, a, a);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorSSEFloat(T (&a)[4]) {
        d = _mm_loadu_ps(&a[0]);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorSSEFloat(T a0, T a1, T a2, T a3) {
        d = _mm_set_ps(a0, a1, a2, a3);
      }

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
        const static __m128 minus_one = _mm_set_ps(-1.0, -1.0, -1.0, -1.0);;
        VectorSSEFloat result;
        result.d = _mm_mul_ps(this->d, minus_one);
        return result;
      }

//      operator float() const {
//        float d0[4];
//        _mm_store_ss(&(d0[0]), d);
//        //alignas(__m128) float d0[4];
//        //_mm_store_ss(&(d0[0]), d);
//        return d0[0];
//      }

      void convert(T(&a)[4]) const {
        _mm_storeu_ps(&a[0], d);
      }
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
  inline VectorSSEFloat fma_plus(VectorSSEFloat a, VectorSSEFloat b, VectorSSEFloat c) {
    VectorSSEFloat d;
    d.d = _mm_fmadd_ps(a.d, b.d, c.d);
    return d;
  }
  inline VectorSSEFloat fma_minus(VectorSSEFloat a, VectorSSEFloat b, VectorSSEFloat c) {
    VectorSSEFloat d;
    d.d = _mm_fmsub_ps(a.d, b.d, c.d);
    return d;
  }
#elif defined(__FMA4__)
  inline VectorSSEFloat fma_plus(VectorSSEFloat a, VectorSSEFloat b, VectorSSEFloat c) {
    VectorSSEFloat d;
    d.d = _mm_macc_ps(a.d, b.d, c.d);
    return d;
  }
  inline VectorSSEFloat fma_minus(VectorSSEFloat a, VectorSSEFloat b, VectorSSEFloat c) {
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
    float a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = std::exp(a_d[i]);
    VectorSSEFloat result(a_d);
#endif
    return result;
  }
  inline VectorSSEFloat sqrt(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
    VectorSSEFloat result;
    result.d = _mm_sqrt_ps(a.d);
#else
    float a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = std::sqrt(a_d[i]);
    VectorSSEFloat result(a_d);
#endif
    return result;
  }
  inline VectorSSEFloat erf(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
    VectorSSEFloat result;
    result.d = _mm_erf_ps(a.d);
#else
    float a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::erf(a_d[i]);
    VectorSSEFloat result(a_d);
#endif
    return result;
  }
  inline VectorSSEFloat erfc(VectorSSEFloat a) {
#if HAVE_INTEL_SVML
    VectorSSEFloat result;
    result.d = _mm_erfc_ps(a.d);
#else
    float a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::erfc(a_d[i]);
    VectorSSEFloat result(a_d);
#endif
    return result;
  }
  //@}

};}; // namespace libint2::simd

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os, libint2::simd::VectorSSEFloat a) {
  float ad[4];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "," << ad[2] << "," << ad[3] << "}";
  return os;
}
//@}

#endif // SSE-only

#ifdef __AVX__

#include <immintrin.h>

namespace libint2 { namespace simd {

  /**
   * SIMD vector of 4 double-precision floating-point real numbers, operations on which use AVX instructions
   * available on recent x86 hardware from Intel (starting with Sandy Bridge processors released in 2011)
   * and AMD (starting with Bulldozer released in 2011).
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
      VectorAVXDouble(T a) {
        d = _mm256_set_pd(a, a, a, a);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorAVXDouble(T (&a)[4]) {
        d = _mm256_loadu_pd(&a[0]);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      VectorAVXDouble(T a0, T a1, T a2, T a3) {
        d = _mm256_set_pd(a0, a1, a2, a3);
      }

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
        const static __m256d minus_one = _mm256_set_pd(-1.0, -1.0, -1.0, -1.0);;
        VectorAVXDouble result;
        result.d = _mm256_mul_pd(this->d, minus_one);
        return result;
      }

//      operator double() const {
//        double d0[4];
//        _mm256_storeu_pd(&(d0[0]), d);
//        // aligned alternative, but efficiency should not matter here
//        //        alignas(__m256d) double d0[4];
//        //        _mm256_store_pd(&(d0[0]), d);
//        return d0[0];
//      }

      void convert(T(&a)[4]) const {
        _mm256_storeu_pd(&a[0], d);
      }
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
  inline VectorAVXDouble fma_plus(VectorAVXDouble a, VectorAVXDouble b, VectorAVXDouble c) {
    VectorAVXDouble d;
    d.d = _mm256_fmadd_pd(a.d, b.d, c.d);
    return d;
  }
  inline VectorAVXDouble fma_minus(VectorAVXDouble a, VectorAVXDouble b, VectorAVXDouble c) {
    VectorAVXDouble d;
    d.d = _mm256_fmsub_pd(a.d, b.d, c.d);
    return d;
  }
#elif defined(__FMA4__)
  inline VectorAVXDouble fma_plus(VectorAVXDouble a, VectorAVXDouble b, VectorAVXDouble c) {
    VectorAVXDouble d;
    d.d = _mm256_facc_pd(a.d, b.d, c.d);
    return d;
  }
  inline VectorAVXDouble fma_minus(VectorAVXDouble a, VectorAVXDouble b, VectorAVXDouble c) {
    VectorAVXDouble d;
    d.d = _mm256_fsub_pd(a.d, b.d, c.d);
    return d;
  }
#endif

  //@}

  //@{ standard functions
  inline VectorAVXDouble exp(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
    VectorAVXDouble result;
    result.d = _mm256_exp_pd(a.d);
#else
    double a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::exp(a_d[i]);
    VectorAVXDouble result(a_d);
#endif
    return result;
  }
  inline VectorAVXDouble sqrt(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
    VectorAVXDouble result;
    result.d = _mm256_sqrt_pd(a.d);
#else
    double a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::sqrt(a_d[i]);
    VectorAVXDouble result(a_d);
#endif
    return result;
  }
  inline VectorAVXDouble erf(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
    VectorAVXDouble result;
    result.d = _mm256_erf_pd(a.d);
#else
    double a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::erf(a_d[i]);
    VectorAVXDouble result(a_d);
#endif
    return result;
  }
  inline VectorAVXDouble erfc(VectorAVXDouble a) {
#if HAVE_INTEL_SVML
    VectorAVXDouble result;
    result.d = _mm256_erfc_pd(a.d);
#else
    double a_d[4]; a.convert(a_d);
    for(int i=0; i<4; ++i) a_d[i] = ::erfc(a_d[i]);
    VectorAVXDouble result(a_d);
#endif
    return result;
  }
  //@}

};}; // namespace libint2::simd

//@{ standard stream operations
inline std::ostream& operator<<(std::ostream& os, libint2::simd::VectorAVXDouble a) {
  double ad[4];
  a.convert(ad);
  os << "{" << ad[0] << "," << ad[1] << "," << ad[2] << "," << ad[3] << "}";
  return os;
}
//@}

#endif // AVX-only

#endif // header guard

