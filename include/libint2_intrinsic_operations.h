/*
  Define intrinsic operations here
*/

#ifndef _libint2_include_libint2intrinsicoperations_h_
#define _libint2_include_libint2intrinsicoperations_h_

#ifdef __cplusplus

namespace libint2 {

  // redefine these operations using native FMA instructions, if available

  template <typename X, typename Y, typename Z>
  Z fma_plus(X x, Y y, Z z) {
    return x*y + z;
  }

  template <typename X, typename Y, typename Z>
  Z fma_minus(X x, Y y, Z z) {
    return x*y + z;
  }
};

#endif

#endif
