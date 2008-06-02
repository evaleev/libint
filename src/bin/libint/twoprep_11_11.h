
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
  /// TwoPRep_11_11_sq is a shell quartet of ERIs
  typedef GenIntegralSet_11_11<CGShell,TwoPRep,mType> TwoPRep_11_11_sq;

  /// TwoPRep_11_11_int is a single ERIs
  typedef GenIntegralSet_11_11<CGF,TwoPRep,mType> TwoPRep_11_11_int;

  /// (ss|ss) shell quartet is not precomputed, but the integral is
  template <>
  inline bool
  GenIntegralSet_11_11<CGF,TwoPRep,mType>::this_precomputed() const
  {
    /// (ss|ss)^{(m)} are precomputed 
#if USE_BRAKET_H
    if (parent_type::bra_.member(0,0).zero() && parent_type::bra_.member(1,0).zero() &&
      parent_type::ket_.member(0,0).zero() && parent_type::ket_.member(1,0).zero())
#else
    if (parent_type::bra_.member(0,0)->zero() && parent_type::bra_.member(1,0)->zero() &&
      parent_type::ket_.member(0,0)->zero() && parent_type::ket_.member(1,0)->zero())
#endif
      return true;
    else
      return false;
  }


};

#endif

