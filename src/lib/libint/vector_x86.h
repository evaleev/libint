
#ifndef _libint2_src_lib_libint_vectorx86_h_
#define _libint2_src_lib_libint_vectorx86_h_

#ifdef __SSE2__

#include <emmintrin.h>

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

      operator double() const {
        double d0;
        _mm_store_sd(&d0, d);
        return d0;
      }

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

};}; // namespace libint2::simd

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

      operator float() const {
        float d0;
        _mm_store_ss(&d0, d);
        return d0;
      }

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

};}; // namespace libint2::simd

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

      operator double() const {
        double d0[4];
        _mm256_store_pd(&(d0[0]), d);
        return d0[0];
      }

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

};}; // namespace libint2::simd

#endif // AVX-only

#endif // header guard

