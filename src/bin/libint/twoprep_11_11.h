/*
 *  Copyright (C) 2004-2017 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_twoprep1111_h_
#define _libint2_src_bin_libint_twoprep1111_h_

#include <integral.h>
#include <integral_11_11.h>

using namespace std;

namespace libint2 {
  
  /**
     Most basic type -- TwoPRep_11_11 --
     has one bfs for each particle in bra and ket.
     Note that GenIntegralSet is initialized with an abstract type libint2::BFSet,
     from which BFS derives.
  */

  /// (ss|ss) shell quartet is not precomputed, but the integral is
  template <>
  inline bool
  GenIntegralSet_11_11<CGF,TwoPRep,mType>::this_precomputed() const
  {
    /// uncontracted (ss|ss)^{(m)} are precomputed
    if (parent_type::bra_.member(0,0).zero() && parent_type::bra_.member(1,0).zero() &&
        parent_type::ket_.member(0,0).zero() && parent_type::ket_.member(1,0).zero() &&
        parent_type::bra_.member(0,0).contracted() == false && parent_type::bra_.member(1,0).contracted() == false &&
        parent_type::ket_.member(0,0).contracted() == false && parent_type::ket_.member(1,0).contracted() == false &&
        parent_type::bra_.member(0,0).deriv().zero() && parent_type::bra_.member(1,0).deriv().zero() &&
        parent_type::ket_.member(0,0).deriv().zero() && parent_type::ket_.member(1,0).deriv().zero()
       )
      return true;
    else
      return false;
  }


};

#endif

