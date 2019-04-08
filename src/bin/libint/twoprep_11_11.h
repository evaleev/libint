/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
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

  namespace detail {
    template <typename Bra>
    bool is_nonderiv_ss_product(Bra&& bra) {
      return bra.member(0,0).zero() && bra.member(1,0).zero() &&
          bra.member(0,0).contracted() == false && bra.member(1,0).contracted() == false &&
          bra.member(0,0).deriv().zero() && bra.member(1,0).deriv().zero();
    };
  }
  /// uncontracted (ss|ss)^{(m)} integral is precomputed (but not shell quartet)
  template <>
  inline bool
  GenIntegralSet_11_11<CGF,TwoPRep,mType>::this_precomputed() const
  {
    if (detail::is_nonderiv_ss_product(parent_type::bra_) && detail::is_nonderiv_ss_product(parent_type::ket_))
      return true;
    else
      return false;
  }

  /// always unroll (ss|ss)^(m) shell set
  template <>
  inline bool
  GenIntegralSet_11_11< CGShell ,TwoPRep,mType>::auto_unroll() const
  {
    if (detail::is_nonderiv_ss_product(parent_type::bra_) && detail::is_nonderiv_ss_product(parent_type::ket_))
      return true;
    else
      return false;
  }

};

#endif

