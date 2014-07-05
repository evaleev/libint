
#ifndef _libint2_src_bin_libint_onep11_h_
#define _libint2_src_bin_libint_onep11_h_

#include <integral.h>
#include <integral_1_1.h>

using namespace std;

namespace libint2 {
  
  /// (s|s) shell quartet is not precomputed, but the integral is
  template <>
  inline bool
  GenIntegralSet_1_1<CGF,OverlapOper,EmptySet>::this_precomputed() const
  {
    /// uncontracted (s|s) are precomputed
    if (parent_type::bra_.member(0,0).zero() &&
        parent_type::ket_.member(0,0).zero() &&
        parent_type::bra_.member(0,0).contracted() == false &&
        parent_type::ket_.member(0,0).contracted() == false &&
        parent_type::bra_.member(0,0).deriv().zero() &&
        parent_type::ket_.member(0,0).deriv().zero()
       )
      return true;
    else
      return false;
  }

  /// (s|T|s) shell quartet is not precomputed, but the integral is
  template <>
  inline bool
  GenIntegralSet_1_1<CGF,KineticOper,EmptySet>::this_precomputed() const
  {
    /// uncontracted (s|s) are precomputed
    if (parent_type::bra_.member(0,0).zero() &&
        parent_type::ket_.member(0,0).zero() &&
        parent_type::bra_.member(0,0).contracted() == false &&
        parent_type::ket_.member(0,0).contracted() == false &&
        parent_type::bra_.member(0,0).deriv().zero() &&
        parent_type::ket_.member(0,0).deriv().zero()
       )
      return true;
    else
      return false;
  }

};

#endif

