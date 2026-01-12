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

#ifndef _libint2_src_lib_libint_vector_h_
#define _libint2_src_lib_libint_vector_h_

#if defined(__cplusplus)

#include <libint2/util/type_traits.h>

#include <algorithm>
#include <cstddef>

namespace libint2 {

/**
   Contains data types that support SIMD-style computation on vectors of
   numbers.
*/
namespace simd {

// add __declspec(align(N*sizeof(T))) ?

/**
 * Vector<N,T> is used by vectorized Libint library as fixed-length vectors
 * amenable for SIMD-style parallelism Vectorization via this class should be
 * the last-resort measure if no specialized implementation is available
 */
template <std::size_t N, typename T>
struct Vector {
  T d[N];

  /**
   * creates a vector of default-initialized values.
   */
  Vector() {}

  /** Initializes all elements to the same value
   *  @param a the value to which all elements will be set
   */
  Vector(T a) { std::fill_n(&(d[0]), N, a); }

  /**
   * creates a vector of values initialized by an ordinary static-sized array
   */
  Vector(T (&a)[N]) { std::copy(&a[0], &a[0] + N, &d[0]); }

  Vector& operator=(T a) {
    for (std::size_t i = 0; i < N; ++i) d[i] = a;
    return *this;
  }

  Vector& operator+=(Vector a) {
    for (std::size_t i = 0; i < N; ++i) d[i] += a.d[i];
    return *this;
  }

  Vector& operator-=(Vector a) {
    for (std::size_t i = 0; i < N; ++i) d[i] -= a.d[i];
    return *this;
  }

  operator double() const { return d[0]; }
};

//@{ arithmetic operators
template <std::size_t N, typename T>
Vector<N, T> operator*(T a, Vector<N, T> b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = a * b.d[i];
  return c;
}

template <std::size_t N, typename T>
Vector<N, T> operator*(Vector<N, T> a, T b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = b * a.d[i];
  return c;
}

template <std::size_t N, typename T>
Vector<N, T> operator*(int a, Vector<N, T> b) {
  if (a == 1)
    return b;
  else {
    Vector<N, T> c;
    for (std::size_t i = 0; i < N; ++i) c.d[i] = a * b.d[i];
    return c;
  }
}

template <std::size_t N, typename T>
Vector<N, T> operator*(Vector<N, T> a, int b) {
  if (b == 1)
    return a;
  else {
    Vector<N, T> c;
    for (std::size_t i = 0; i < N; ++i) c.d[i] = b * a.d[i];
    return c;
  }
}

template <std::size_t N, typename T>
Vector<N, T> operator*(Vector<N, T> a, Vector<N, T> b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = a.d[i] * b.d[i];
  return c;
}

template <std::size_t N, typename T>
Vector<N, T> operator+(Vector<N, T> a, Vector<N, T> b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = a.d[i] + b.d[i];
  return c;
}

template <std::size_t N, typename T>
Vector<N, T> operator-(Vector<N, T> a, Vector<N, T> b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = a.d[i] - b.d[i];
  return c;
}

template <std::size_t N, typename T>
Vector<N, T> operator/(Vector<N, T> a, Vector<N, T> b) {
  Vector<N, T> c;
  for (std::size_t i = 0; i < N; ++i) c.d[i] = a.d[i] / b.d[i];
  return c;
}

//@}

};  // namespace simd
};  // namespace libint2

namespace libint2 {

template <std::size_t N, typename T>
struct is_vector<simd::Vector<N, T> > {
  static const bool value = true;
};

template <std::size_t N, typename T>
struct vector_traits<simd::Vector<N, T> > {
  typedef T value_type;
  static const std::size_t extent = N;
};

}  // namespace libint2

#include "vector_ppc.h"
#include "vector_x86.h"

#endif  // C++ only

#endif
