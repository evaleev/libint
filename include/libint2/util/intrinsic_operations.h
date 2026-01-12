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

#ifndef _libint2_include_libint2intrinsicoperations_h_
#define _libint2_include_libint2intrinsicoperations_h_

#include <libint2/config.h>

#ifdef __cplusplus

namespace libint2 {

//@{ Floating-point-Multiply-Add (FMA) instructions. Redefine these operations
// using native FMA instructions, if available (see, e.g. vector_x86.h)

#if defined(LIBINT_GENERATE_FMA)
#if defined(LIBINT_HAS_CXX11)
/// @return x*y+z
template <typename X, typename Y, typename Z>
inline auto fma_plus(X x, Y y, Z z) -> decltype(x * y + z) {
  return x * y + z;
}

/// @return x*y-z
template <typename X, typename Y, typename Z>
inline auto fma_minus(X x, Y y, Z z) -> decltype(x * y - z) {
  return x * y - z;
}
#else  // LIBINT_HAS_CXX11
#error "support for FMA requires compiler capable of C++11 or later"
#endif  // LIBINT_HAS_CXX11
#endif  // LIBINT_GENERATE_FMA

//@}

};  // namespace libint2

/**
   these macros define bzero, copy, and inc operations:
*/
/** X[i] = 0 */
#define _libint2_static_api_bzero_short_(X, nelem) \
  for (int i = 0; i < (nelem); ++i) {              \
    (X)[i] = 0.0;                                  \
  }
/** X[i] = Y[i] */
#define _libint2_static_api_copy_short_(X, Y, nelem) \
  for (int i = 0; i < (nelem); ++i) {                \
    (X)[i] = (Y)[i];                                 \
  }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_short_(X, Y, nelem, a) \
  for (int i = 0; i < (nelem); ++i) {                    \
    (X)[i] = (a) * (Y)[i];                               \
  }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_vec_short_(X, Y, nelem, a, vl) \
  for (int i = 0, iv = 0; i < (nelem) / (vl); ++i) {             \
    for (int v = 0; v < (vl); ++v, ++iv) {                       \
      (X)[iv] = (a)[v] * (Y)[iv];                                \
    }                                                            \
  }
/** X[i] += a*Y[i] */
#if defined(LIBINT_GENERATE_FMA)
#define _libint2_static_api_inc_short_(X, Y, nelem, a) \
  for (int i = 0; i < (nelem); ++i) {                  \
    (X)[i] = libint2::fma_plus((a), (Y)[i], (X)[i]);   \
  }
#else
#define _libint2_static_api_inc_short_(X, Y, nelem, a) \
  for (int i = 0; i < (nelem); ++i) {                  \
    (X)[i] += (a) * (Y)[i];                            \
  }
#endif
/** X[i] += Y[i] */
#define _libint2_static_api_inc1_short_(X, Y, nelem) \
  for (int i = 0; i < (nelem); ++i) {                \
    (X)[i] += (Y)[i];                                \
  }

#endif

#endif
