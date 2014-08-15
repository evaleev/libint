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

#ifndef _libint2_include_libint2intrinsicoperations_h_
#define _libint2_include_libint2intrinsicoperations_h_

#ifdef __cplusplus

namespace libint2 {

  //@{ Floating-point-Multiply-Add (FMA) instructions. Redefine these operations using native FMA instructions, if available (see, e.g. vector_x86.h)

  /// @return x*y+z
  template <typename X, typename Y, typename Z>
  Z fma_plus(X x, Y y, Z z) {
    return x*y + z;
  }

  /// @return x*y-z
  template <typename X, typename Y, typename Z>
  Z fma_minus(X x, Y y, Z z) {
    return x*y - z;
  }

  //@}

};

/**
   these macros define bzero, copy, and inc operations:
*/
/** X[i] = 0 */
#define _libint2_static_api_bzero_short_(X,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] = 0.0; }
/** X[i] = Y[i] */
#define _libint2_static_api_copy_short_(X,Y,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] = (Y)[i]; }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_short_(X,Y,nelem,a) for(int i=0; i < (nelem); ++i) { (X)[i] = (a) * (Y)[i]; }
/** X[i] = a*Y[i] */
#define _libint2_static_api_scale_vec_short_(X,Y,nelem,a,vl) for(int i=0, iv=0; i < (nelem)/(vl); ++i) { for(int v=0; v < (vl); ++v, ++iv) { (X)[iv] = (a)[v] * (Y)[iv]; } }
/** X[i] += a*Y[i] */
#define _libint2_static_api_inc_short_(X,Y,nelem,a) for(int i=0; i < (nelem); ++i) { (X)[i] = libint2::fma_plus((a),(Y)[i],(X)[i]); }
/** X[i] += Y[i] */
#define _libint2_static_api_inc1_short_(X,Y,nelem) for(int i=0; i < (nelem); ++i) { (X)[i] += (Y)[i]; }

#endif

#endif
