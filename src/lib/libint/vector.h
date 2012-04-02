
#ifndef _libint2_src_lib_libint_vector_h_
#define _libint2_src_lib_libint_vector_h_

#include <unistd.h>

#if defined(__cplusplus)

#include <algorithm>

namespace libint2 {

  /**
     Contains data types that support SIMD-style computation on vectors of numbers.
  */
  namespace simd {

  // add __declspec(align(N*sizeof(T))) ?

  /**
   * Vector<N,T> is used by vectorized Libint library as fixed-length vectors amenable for SIMD-style parallelism
   * Vectorization via this class should be the last-resort measure if no specialized implementation is available
   */
  template <size_t N, typename T>
    struct Vector {

      T d[N];

      /**
       * creates a vector of default-initialized values.
       */
      Vector() {}

      /** Initializes all elements to the same value
       *  @param a the value to which all elements will be set
       */
      Vector(T a) {
        std::fill_n(&(d[0]), N, a);
      }

      /**
       * creates a vector of values initialized by an ordinary static-sized array
       */
      Vector(T (&a)[N]) {
        std::copy(&a[0], &a[0]+N, &d[0]);
      }

      Vector& operator=(T a) {
        for(size_t i=0; i<N; ++i)
          d[i] = a;
        return *this;
      }

      Vector& operator+=(Vector a) {
        for(size_t i=0; i<N; ++i)
          d[i] += a.d[i];
        return *this;
      }

      Vector& operator-=(Vector a) {
        for(size_t i=0; i<N; ++i)
          d[i] -= a.d[i];
        return *this;
      }

      operator double() const {
        return d[0];
      }

  };

  //@{ arithmetic operators
  template <size_t N, typename T>
  Vector<N,T> operator*(T a, Vector<N,T> b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = a * b.d[i];
    return c;
  }

  template <size_t N, typename T>
  Vector<N,T> operator*(Vector<N,T> a, T b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = b * a.d[i];
    return c;
  }

  template <size_t N, typename T>
  Vector<N,T> operator*(int a, Vector<N,T> b) {
    if (a == 1)
      return b;
    else {
      Vector<N, T> c;
      for (size_t i = 0; i < N; ++i)
        c.d[i] = a * b.d[i];
      return c;
    }
  }

  template <size_t N, typename T>
  Vector<N,T> operator*(Vector<N,T> a, int b) {
    if (b == 1)
      return a;
    else {
      Vector<N, T> c;
      for (size_t i = 0; i < N; ++i)
        c.d[i] = b * a.d[i];
      return c;
    }
  }

  template <size_t N, typename T>
  Vector<N,T> operator*(Vector<N,T> a, Vector<N,T> b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = a.d[i] * b.d[i];
    return c;
  }

  template <size_t N, typename T>
  Vector<N,T> operator+(Vector<N,T> a, Vector<N,T> b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = a.d[i] + b.d[i];
    return c;
  }

  template <size_t N, typename T>
  Vector<N,T> operator-(Vector<N,T> a, Vector<N,T> b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = a.d[i] - b.d[i];
    return c;
  }

  template <size_t N, typename T>
  Vector<N,T> operator/(Vector<N,T> a, Vector<N,T> b) {
    Vector<N,T> c;
    for(size_t i=0; i<N; ++i)
      c.d[i] = a.d[i] / b.d[i];
    return c;
  }


  //@}

};}; // namespace libint2::simd

#include <vector_x86.h>
#include <vector_ppc.h>

#endif // C++ only

#endif
