
#ifndef _libint2_src_lib_libint_vector_h_
#define _libint2_src_lib_libint_vector_h_

#include <unistd.h>

namespace libint2 {

  // add __declspec(align(N*sizeof(T))) ?
  template <size_t N, typename T>
    struct Vector {

      T d[N];

      Vector() {}

      Vector(T a) {
        for(size_t i=0; i<N; ++i)
          d[i] = a;
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

};

# include <vector_x86.h>

#endif
